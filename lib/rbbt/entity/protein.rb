require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/statistics/hypergeometric'
require 'rbbt/network/paths'
require 'rbbt/entity/gene'

Workflow.require_workflow "Translation"

module Protein
  extend Entity
  include Entity::Adjacent
  include Entity::Enriched

  self.annotation :format
  self.annotation :organism

  self.format = "Ensembl Protein ID"

  def self.ensp2sequence(organism, protein)
    @ensp2sequence ||= {}
    @ensp2sequence[organism] ||= Organism.protein_sequence(organism).tsv :persist => true
    if Array === protein
      @ensp2sequence[organism].values_at *protein
    else
      @ensp2sequence[organism][protein]
    end
  end

  def self.ensp2enst(organism, protein)
    @ensp2enst ||= {}
    @ensp2enst[organism] ||= Organism.transcripts(organism).tsv(:type => :single, :key_field => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :persist => true)
    @ensp2enst[organism][protein]
  end

  property :uniprot => :array2single do
    @uniprot ||= to "UniProt/SwissProt Accession"
  end

  property :ensembl => :array2single do
    @ensembl = to "Ensembl Protein ID"
  end

  property :transcript => :array2single do
    ensembl.collect{|ensp|
      Protein.ensp2enst(organism, ensp)
    }
  end

  property :ensembl_protein_image_url => :single2array do
    ensembl_url = if organism == "Hsa" then "www.ensembl.org" else "#{organism.sub(/.*\//,'')}.archive.ensembl.org" end
    "http://#{ensembl_url}/Homo_sapiens/Component/Transcript/Web/TranslationImage?db=core;p=#{ensembl};_rmd=d2a8;export=svg"
  end

  property :to! => :array2single do |new_format|
    return self if format == new_format
    Protein.setup(Translation.job(:translate_protein, "", :organism => organism, :proteins => self, :format => new_format).exec, new_format, organism)
  end

  property :to => :array2single do |new_format|
    return self if format == new_format
    to!(new_format).collect!{|v| v.nil? ? nil : v.first}
  end

  property :gene do
    Gene.setup(to("Ensembl Protein ID").clean_annotations, "Ensembl Protein ID", organism)
  end

  property :pfam => :array2single do
    index = Organism.gene_pfam(organism).tsv :flat, :persist => true
    index.unnamed = true
    pfam = index.values_at(*self).flatten
    Pfam.setup pfam
  end

  property :sequence => :array2single do
    Protein.ensp2sequence(organism, self.ensembl)
  end

  property :sequence_length => :array2single do
    sequence.collect{|seq| seq.nil? ? nil : seq.length}
  end
end

