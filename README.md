# deeplearningmetagenomics
Understanding relationships between complex combinations of OTUs and chronic disease


## Day 1 feedback
- why use OTUs?
- we can instead just use 16S reads to map to whole genomes
- can feed the full gff3 file to PICRUST
- use HMMer to convert whole genomes to sequences of domains
- essentially like a corpus where each domain = a word
- feed these to a transformer model to train
  - we can train the encoder
  - if we have matched input:output, can separately train the decoder side