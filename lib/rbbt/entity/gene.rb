require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/sources/entrez'

Workflow.require_workflow "Translation"

module Gene
  extend Entity

  self.annotation :format
  self.annotation :organism

  self.format = Organism::Hsa.identifiers.all_fields

  def name
    to "Associated Gene Name"
  end

  def long_name
    if Array === self
      to("Entrez Gene ID").collect{|id| gene = Entrez.get_gene(id); gene.nil? ? nil : gene.description}
    else
      gene = Entrez.get_gene(to("Entrez Gene ID"))
      gene.nil? ? nil : gene.description
    end
  end


  def description
    if Array === self
      to("Entrez Gene ID").collect{|id| gene = Entrez.get_gene(id); gene.nil? ? nil : gene.summary}
    else
      gene = Entrez.get_gene(to("Entrez Gene ID"))
      gene.nil? ? nil : gene.summary
    end
  end

  def to!(new_format)
    if Array === self
      Gene.setup(Translation.job(:tsv_translate, "", :organism => organism, :genes => self, :format => new_format).exec.values_at(*self), new_format, organism)
    else
      Gene.setup(Translation.job(:tsv_translate, "", :organism => organism, :genes => [self], :format => new_format).exec[self], new_format, organism)
    end
  end

  def to(new_format)
    return self if format == new_format
    if Array === self
      to!(new_format).collect!{|v| v.nil? ? nil : v.first}
    else
      v = to!(new_format)
      v.nil? ? nil : v.first
    end
  end

  def self2pfam
    index = Organism.gene_pfam(organism).tsv :type => :flat, :persist => true
    if Array === self
      index.values_at(*self).flatten
    else
      index[self]
    end
  end

  def chromosome
    chr = Organism.gene_positions(organism).tsv :fields => ["Chromosome Name"], :type => :single, :persist => true
    if Array === self
      to("Ensembl Gene ID").collect do |gene|
        chr[gene]
      end
    else
      chr[to("Ensembl Gene ID")]
    end
  end

  def range
    pos = Organism.gene_positions(organism).tsv :fields => ["Gene Start", "Gene End"], :type => :list, :persist => true, :cast => :to_i
    if Array === self
      to("Ensembl Gene ID").collect do |gene|
        next if not pos.include? gene
        Range.new *pos[gene]
      end
    else
      return nil if not pos.include? to("Ensembl Gene ID")
      Range.new *pos[to("Ensembl Gene ID")]
    end
  end

  def self2transcripts
    if Array === self
      gene_transcripts = Organism.gene_transcripts(organism).tsv :persit => true
      gene_transcripts.select self.to "Ensembl Gene ID"
    else
      self.make_list.self2transcripts[self.to "Ensembl Gene ID"]
    end
  end

  def self2max_protein_length
    if Array === self
      self2transcripts = self.self2transcripts
      protein_sequence = Organism.protein_sequence(organism).tsv :persist => true
      transcript_proteins = Organism.transcripts(organism).tsv :fields => ["Ensembl Protein ID"], :type => :single, :persist => true
      
      Misc.process_to_hash(self.to "Ensembl Gene ID") do |list|
        self2transcripts.values_at(*list).collect{|trans| 
          proteins = transcript_proteins.values_at(*trans).compact
          protein_sequence.values_at(*proteins).collect{|seq| 
           seq.nil? ? 0 : seq.length
          }.max
        }
      end
    else
      self.make_list.self2max_protein_length[self.to "Ensembl Gene ID"]
    end
  end

  def self2max_transcript_length
    if Array === self
      self2transcripts = self.self2transcripts
      transcript_sequence = Organism.transcript_sequence(organism).tsv
      
      Misc.process_to_hash(self.to "Ensembl Gene ID") do |list|
        self2transcripts.values_at(*list).collect{|trans| 
          transcript_sequence.values_at(*trans).collect{|seq| 
           seq.nil? ? 0 : seq.length
          }.max
        }
      end
    else
      self.make_list.self2max_transcript_length[self.to "Ensembl Gene ID"]
    end
  end

end

module Transcript
  extend Entity

  def to!(new_format)
    if Array === self
      Gene.setup(Translation.job(:tsv_probe_translate, "", :organism => organism, :genes => self, :format => new_format).exec.values_at(*self), new_format, organism)
    else
      Gene.setup(Translation.job(:tsv_probe_translate, "", :organism => organism, :genes => [self], :format => new_format).exec[self], new_format, organism)
    end
  end

  def to(new_format)
    return self if format == new_format
    if Array === self
      to!(new_format).collect{|v| v.nil? ? nil : v.first}
    else
      v = to!(new_format)
      v.nil? ? nil : v.first
    end
  end
end

