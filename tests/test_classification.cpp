#define main uniann_program_main
#include "../src/uniann.cpp"
#undef main

#include <cassert>

static TranscriptAnnotation transcript(
    const vector<GenomicInterval> &exons,
    const vector<SpliceJunction> &junctions
) {
    TranscriptAnnotation result;
    result.sequence_id = "synthetic";
    result.genomic_start_0based = exons.front().start_0based;
    result.genomic_end_0based_exclusive = exons.back().end_0based_exclusive;
    result.exons = exons;
    result.cds_intervals = exons;
    result.junctions = junctions;
    return result;
}

int main() {
    const auto reference = transcript(
        {{0, 50}, {100, 150}, {200, 250}},
        {{50, 98, '+'}, {150, 198, '+'}});

    const auto alternative_donor = transcript(
        {{0, 55}, {100, 150}, {200, 250}},
        {{55, 98, '+'}, {150, 198, '+'}});
    assert(classify_splice_difference(reference, alternative_donor) ==
           "alternative_donor");

    const auto alternative_acceptor = transcript(
        {{0, 50}, {109, 150}, {200, 250}},
        {{50, 107, '+'}, {150, 198, '+'}});
    assert(classify_splice_difference(reference, alternative_acceptor) ==
           "alternative_acceptor");

    const auto skipped = transcript(
        {{0, 50}, {200, 250}}, {{50, 198, '+'}});
    assert(classify_splice_difference(reference, skipped) == "exon_skipping");
    assert(classify_splice_difference(skipped, reference) == "alternative_exon");

    const auto retained = transcript(
        {{0, 150}, {200, 250}}, {{150, 198, '+'}});
    assert(classify_splice_difference(reference, retained) ==
           "intron_retention");

    const auto complex = transcript(
        {{3, 48}, {111, 151}, {205, 248}},
        {{48, 109, '+'}, {151, 203, '+'}});
    assert(classify_splice_difference(reference, complex) == "complex");

    vector<int> positions(103);
    vector<int> states(103, 0);
    for (int i = 0; i < 103; ++i) positions[i] = i;
    for (int i = 2; i <= 40; ++i) states[i] = 1;
    for (int i = 41; i <= 70; ++i) states[i] = 4;
    for (int i = 71; i <= 101; ++i) states[i] = 3;
    const auto converted = build_transcript_annotations(
        positions, states, "coordinate_test");
    assert(converted.size() == 1);
    assert((converted[0].exons == vector<GenomicInterval>{{0, 40}, {72, 103}}));
    assert(converted[0].junctions.size() == 1);
    assert(converted[0].junctions[0].donor_0based == 40);
    assert(converted[0].junctions[0].acceptor_0based == 70);

    const auto phase_test = transcript(
        {{0, 50}, {100, 160}, {200, 260}},
        {{50, 98, '+'}, {160, 198, '+'}});
    ostringstream gff;
    write_hierarchical_transcript(gff, phase_test, "g1", "g1.t1");
    assert(gff.str().find("\t1\tID=g1.t1.cds2") != string::npos);
    assert(gff.str().find("\t1\tID=g1.t1.cds3") != string::npos);

    const int length = 180;
    vector<char> sequence(length, 'A');
    vector<array<double, NUM_STATES>> emissions(length);
    for (auto &row : emissions) row.fill(0.0);
    vector<double> gt(length, NEG_INF), ag(length, NEG_INF);
    vector<double> atg(length, NEG_INF), stop(length, NEG_INF);
    sequence[30] = 'A'; sequence[31] = 'T'; sequence[32] = 'G';
    sequence[60] = 'G'; sequence[61] = 'T';
    sequence[90] = 'A'; sequence[91] = 'G';
    sequence[120] = 'T'; sequence[121] = 'A'; sequence[122] = 'A';
    atg[30] = 5.0; gt[60] = 7.0; ag[90] = 11.0; stop[120] = 13.0;
    const auto transitions = init_transitions();
    const ModelInputs inputs{emissions, gt, ag, atg, stop, sequence, transitions};

    PathMetadata intergenic;
    intergenic.inter_len = 30;
    intergenic.predecessor_state = 0;
    const auto valid_start = evaluate_transition(32, 0, 1, intergenic, inputs);
    assert(valid_start.allowed && valid_start.transition_score == 5.0);
    assert(!evaluate_transition(32, 0, 2, intergenic, inputs).allowed);
    sequence[30] = 'C';
    assert(!evaluate_transition(32, 0, 1, intergenic, inputs).allowed);
    sequence[30] = 'A';
    PathMetadata short_intergenic = intergenic;
    short_intergenic.inter_len = 1;
    assert(!evaluate_transition(32, 0, 1, short_intergenic, inputs).allowed);

    PathMetadata exon;
    exon.exon_len = 40;
    exon.exon_from = 4;
    const auto valid_donor = evaluate_transition(61, 1, 4, exon, inputs);
    assert(valid_donor.allowed && valid_donor.transition_score == 7.0);
    assert(!evaluate_transition(61, 1, 5, exon, inputs).allowed);
    PathMetadata short_exon = exon;
    short_exon.exon_len = 10;
    assert(!evaluate_transition(61, 1, 4, short_exon, inputs).allowed);

    PathMetadata intron;
    intron.intron_len = 28; // plus the AG gives a 30-base intron
    const auto valid_acceptor = evaluate_transition(91, 4, 1, intron, inputs);
    assert(valid_acceptor.allowed && valid_acceptor.transition_score == 11.0);
    assert(!evaluate_transition(91, 4, 2, intron, inputs).allowed);
    PathMetadata short_intron = intron;
    short_intron.intron_len = 20;
    assert(!evaluate_transition(91, 4, 1, short_intron, inputs).allowed);

    PathMetadata terminal_exon;
    terminal_exon.exon_len = 120;
    terminal_exon.exon_from = 4;
    const auto valid_stop = evaluate_transition(122, 1, 0, terminal_exon, inputs);
    assert(valid_stop.allowed && valid_stop.transition_score == 13.0);
    assert(!evaluate_transition(122, 2, 0, terminal_exon, inputs).allowed);
    sequence[120] = 'C';
    assert(!evaluate_transition(122, 1, 0, terminal_exon, inputs).allowed);

    string geometry_reason;
    const auto too_short_single = transcript({{0, 90}}, {});
    assert(!valid_transcript_geometry(too_short_single, geometry_reason));
    assert(geometry_reason == "single_exon_minimum_length_violation");

    return 0;
}
