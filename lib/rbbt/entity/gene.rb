require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/sources/entrez'
require 'rbbt/entity/protein'
require 'rbbt/entity/pmid'

Workflow.require_workflow "Translation"

module Gene
  extend Entity

  self.annotation :format
  self.annotation :organism

  self.format = Organism::Hsa.identifiers.all_fields - ["Ensembl Protein ID", "Ensembl Transcript ID"]

  property :to! => :array2single do |new_format|
    return self if format == new_format
    Gene.setup(Translation.job(:tsv_translate, "", :organism => organism, :genes => self, :format => new_format).exec.values_at(*self), new_format, organism)
  end

  property :to => :array2single do |new_format|
    return self if format == new_format
    to!(new_format).collect!{|v| v.nil? ? nil : v.first}
  end

  property :ensembl => :array2single do
    @ensembl ||= to "Ensembl Gene ID"
  end

  property :entrez => :array2single do
    @entrez ||= to "Entrez Gene ID"
  end


  property :name => :array2single do
    @name ||= to "Associated Gene Name"
  end

  property :chr_start => :array2single do
    @chr_start = begin
                   Organism.gene_positions(organism).tsv(:persist => true, :type => :single, :cast => :to_i, :fields => ["Gene Start"]).values_at *self
                 end
  end

  property :go_bp_terms => :array2single do
    @go_bp_terms ||= Organism.gene_go_bp(organism).tsv(:persist => true, :key_field => "Ensembl Gene ID", :fields => ["GO ID"], :type => :flat).values_at *self.ensembl
  end

  property :long_name => :single2array do
    gene = Entrez.get_gene(to("Entrez Gene ID"))
    gene.nil? ? nil : gene.description.flatten.first
  end

  property :description => :single2array do
    gene = Entrez.get_gene(to("Entrez Gene ID"))
    gene.nil? ? nil : gene.summary.flatten.first
  end

  property :transcripts => :array2single do
    gene_transcripts = Organism.gene_transcripts(organism).tsv :persist => true
    gene_transcripts.unnamed = true
    res = gene_transcripts.values_at(*self.ensembl)
    res.each{|l| Transcript.setup(l, "Ensembl Transcript ID", organism)}
    res
  end

  property :proteins  => :array2single do
    @proteins ||= begin
                    transcripts = self.transcripts
                    all_transcripts = Transcript.setup(transcripts.flatten, "Ensembl Transcript ID", organism)
                    transcript2protein = nil

                    transcript2protein = Misc.process_to_hash(all_transcripts){|list|
                      list.protein
                    }

                    res = nil
                    res = transcripts.collect{|list|
                      Protein.setup(transcript2protein.values_at(*list), "Ensembl Protein ID", organism)
                    }

                    res.each{|l| 
                    }
                    res
                  end
  end

  property :max_transcript_length => :array2single do
    transcripts.collect{|list| list.sequence_length.compact.max}
  end

  property :max_protein_length => :array2single do
    @max_protein_length ||= begin
                              proteins = self.proteins
                              all_proteins = Protein.setup(proteins.flatten, "Ensembl Protein ID", organism)
                              lengths = Misc.process_to_hash(all_proteins){|list| list.sequence_length}
                              proteins.collect{|list| lengths.values_at(*list).compact.max}
                            end
  end

  property :chromosome => :array2single do
    chr = Organism.gene_positions(organism).tsv :fields => ["Chromosome Name"], :type => :single, :persist => true
    chr.unnamed = true
    if Array === self
      to("Ensembl Gene ID").collect do |gene|
        chr[gene]
      end
    else
      chr[to("Ensembl Gene ID")]
    end
  end

  property :range => :array2single do
    pos = Organism.gene_positions(organism).tsv :fields => ["Gene Start", "Gene End"], :type => :list, :persist => true, :cast => :to_i
    to("Ensembl Gene ID").collect do |gene|
      next if not pos.include? gene
      Range.new *pos[gene]
    end
  end

  property :articles => :array2single do
    @articles ||= begin
                    PMID.setup(Organism.gene_pmids(organism).tsv(:persist => true, :fields => ["PMID"], :type => :flat).values_at *self.entrez)
                  end 
  end
end

module Transcript
  extend Entity

  self.annotation :format
  self.annotation :organism

  self.format = "Ensembl Transcript ID"

  property :to! => :array2single do |new_format|
    return self if format == new_format
    Gene.setup(Translation.job(:tsv_probe_translate, "", :organism => organism, :genes => self, :format => new_format).exec.values_at(*self), new_format, organism)
  end

  property :to => :array2single do |new_format|
    return self if format == new_format
    to!(new_format).collect!{|v| v.nil? ? nil : v.first}
  end

  def ensembl
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

  property :protein  => :array2single do
    transcript_protein = Organism.transcripts(organism).tsv :single, :persist => true, :fields => ["Ensembl Protein ID"]
    transcript_protein.unnamed = true

    res = transcript_protein.values_at(*self.ensembl)
    Protein.setup(res, "Ensembl Protein ID", organism)
    res
  end
end

