// Cohort stage: build a normal panel per caller, filter each pair's calls, and merge the callers.
// Input rows are (id, normal_id, caller, role, vcf). Callers are whatever values appear in the
// caller column; nothing here is specific to a particular caller.
//   role tumor  - a pair's tumour calls from one caller
//   role normal - the pair's normal calls from one caller, used for the panel and for subtraction
//   role panel  - a frozen panel for one caller, used instead of building one

include { PANEL_MERGE; PANEL_FILTER } from './jasmine_panel.nf'
include { TUMOUR_FILTER; NORMAL_ALT_ONLY } from './tumour_filter.nf'
include { NORMALIZE_VCF; TUMOUR_CONSENSUS_VCF; TUMOUR_CONSENSUS_CALL } from './tumour_consensus.nf'

workflow SV_COHORT {

    take:
    rows

    main:
    def no_vcf = file("${projectDir}/assets/NO_FILE")
    def min_support = params.panel_min_support as int
    def use_matched = params.matched_normal_subtract.toString().toBoolean()

    tumours = rows.filter { it[3] == 'tumor' }.map { id, nid, caller, role, vcf -> tuple(id, caller, vcf) }
    normal_rows = rows.filter { it[3] == 'normal' }.map { id, nid, caller, role, vcf -> tuple(id, nid, caller, vcf) }
    // usable is false when a normal VCF has no way to tell whether the normal carries each SV
    def normals = normal_rows.map { id, nid, caller, vcf -> tuple(id, nid, caller, vcf, true) }
    if (params.normal_alt_only.toString().toBoolean()) {
        def expressions = params.normal_filters ?: [:]
        NORMAL_ALT_ONLY(normal_rows.map { id, nid, caller, vcf -> tuple(id, nid, caller, vcf, (expressions[caller] ?: '').toString()) })
        normals = NORMAL_ALT_ONLY.out.vcf.map { id, nid, caller, vcf, flag -> tuple(id, nid, caller, vcf, flag.text.trim() == 'true') }
    }
    frozen = rows.filter { it[3] == 'panel' }.map { id, nid, caller, role, vcf -> tuple(caller, vcf) }

    // One VCF per distinct normal and caller. A normal shared by several pairs is counted once,
    // using the VCF from the first of those pairs by id.
    distinct_normals = normals
        .map { id, nid, caller, vcf, usable -> tuple([caller, nid], id, vcf) }
        .groupTuple()
        .map { key, ids, vcfs ->
            def first = ids.indexOf(ids.min())
            tuple(key[0], vcfs[first])
        }
        .groupTuple()
        .map { caller, vcfs -> tuple(caller, vcfs.sort { it.name }) }

    // Callers with a frozen panel skip the build
    // toList().map { [it] } keeps the list as a single item when combined
    frozen_callers = frozen.map { it[0] }.toList().map { [it] }
    to_build = distinct_normals
        .combine(frozen_callers)
        .filter { caller, vcfs, fixed -> !fixed.contains(caller) }
        .map { caller, vcfs, fixed -> tuple(caller, vcfs, min_support) }

    PANEL_MERGE(to_build)
    PANEL_FILTER(PANEL_MERGE.out.merged)
    panels = PANEL_FILTER.out.panel.mix(frozen).toList().map { [it] }

    // Pair each tumour VCF with its caller's panel and the pair's own normal VCF
    tumour_keyed = tumours
        .combine(panels)
        .map { id, caller, vcf, panel_list ->
            def panel = panel_list.find { it[0] == caller }
            tuple([id, caller], 'T', vcf, panel ? panel[1] : no_vcf, panel ? true : false)
        }
    normal_keyed = normals
        .map { id, nid, caller, vcf, usable -> tuple([id, caller], 'N', vcf, no_vcf, usable) }

    to_filter = tumour_keyed
        .mix(normal_keyed)
        .groupTuple()
        .filter { key, kinds, vcfs, panel_vcfs, panel_flags -> kinds.contains('T') }
        .map { key, kinds, vcfs, panel_vcfs, panel_flags ->
            def t = kinds.indexOf('T')
            def n = kinds.indexOf('N')
            def matched = use_matched && n >= 0 && panel_flags[n]
            tuple(key[0], key[1], vcfs[t], panel_vcfs[t], panel_flags[t], matched ? vcfs[n] : no_vcf, matched)
        }

    TUMOUR_FILTER(to_filter)
    NORMALIZE_VCF(TUMOUR_FILTER.out.filtered_vcf)

    // One merge per pair, with the callers in a fixed order
    per_pair = NORMALIZE_VCF.out.vcf
        .groupTuple()
        .map { id, callers, vcfs ->
            def order = (0..<callers.size()).collect { it }.sort { callers[it] }
            tuple(id, order.collect { callers[it] }, order.collect { vcfs[it] })
        }

    TUMOUR_CONSENSUS_VCF(per_pair)
    TUMOUR_CONSENSUS_CALL(TUMOUR_CONSENSUS_VCF.out.id, TUMOUR_CONSENSUS_VCF.out.consensus_vcf, TUMOUR_CONSENSUS_VCF.out.n_callers)
}
