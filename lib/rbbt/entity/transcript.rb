require 'rbbt/entity'
require 'rbbt/entity/gene'

module Transcript
  extend Entity

  self.annotation :format
  self.annotation :organism

  self.format = "Ensembl Transcript ID"

  def self.enst2ensg(organism, transcript)
    @@enst2ensg ||= {}
    @@enst2ensg[organism] ||= Organism.gene_transcripts(organism).tsv(:type => :single, :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Gene ID"], :persist => true).tap{|o| o.unnamed = true}
    res = if Array === transcript
            @@enst2ensg[organism].values_at *transcript
          else
            @@enst2ensg[organism][transcript]
          end
    Gene.setup(res, "Ensembl Gene ID", organism)
  end

  def self.enst2ensp(organism, transcript)
    @@enst2ensp ||= {}
    @@enst2ensp[organism] ||= Organism.transcripts(organism).tsv(:type => :single, :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Protein ID"], :persist => true)
    res = if Array === transcript
            @@enst2ensp[organism].values_at *transcript
          else
            @@enst2ensp[organism][transcript]
          end
    Protein.setup(res, "Ensembl Protein ID", organism)
  end


  property :to! => :array2single do |new_format|
    return self if format == new_format
    Gene.setup(Translation.job(:tsv_translate_probe, "", :organism => organism, :probes => self, :format => new_format).exec.values_at(*self), new_format, organism)
  end

  property :to => :array2single do |new_format|
    return self if format == new_format
    to!(new_format).collect!{|v| v.nil? ? nil : v.first}
  end

  property :ensembl => :array2single do
    to "Ensembl Transcript ID"
  end

  property :sequence => :array2single do
    transcript_sequence = Organism.transcript_sequence(organism).tsv :persist => true
    transcript_sequence.unnamed = true
    transcript_sequence.values_at *self.ensembl
  end

  property :sequence_length => :array2single do
    sequence.collect{|s|
      s.nil? ? nil : s.length
    }
  end

  property :gene => :array2single do
    Transcript.enst2ensg(organism, self)
  end

  property :protein => :array2single do
    Transcript.enst2ensp(organism, self)
  end
end
