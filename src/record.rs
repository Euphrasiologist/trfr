/// The record fields in a `trf` table
#[derive(Debug, Default, PartialEq)]
pub struct Record {
    /// The name of the fasta record
    pub seq_id: String,
    /// Start index of the repeat
    pub start: usize,
    /// End index of the repeat
    pub end: usize,
    /// The period of the repeat
    pub period: u16,
    /// Number of copies aligned with the consensus pattern
    pub copy_number: f32,
    /// Size of consensus pattern (may differ slightly from the period size).
    pub consensus_pattern_size: u16,
    /// Percent of matches between adjacent copies overall
    pub perc_matches: u8,
    /// Percent of indels between adjacent copies overall
    pub perc_indels: u8,
    /// Alignment score
    pub alignment_score: u32,
    /// Percentages of the nucleotides
    pub perc_a: u8,
    pub perc_c: u8,
    pub perc_g: u8,
    pub perc_t: u8,
    /// Entropy measure based on percent composition
    pub entropy: f32,
    /// The repeat pattern itself
    pub consensus_pattern: String,
    /// The longer repeat sequence extracted from the reads
    pub repeat_seq: String,
}
