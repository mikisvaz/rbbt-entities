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

  def ensembl
    to "Ensembl Protein ID"
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
    Organism.gene_pfam(organism).tsv :flat, :persist => true
    pfam = index.values_at(*self).flatten
    Pfam.setup pfam
  end

  property :sequence => :array2single do
    @protein_sequence ||= begin
                            protein_sequence = Organism.protein_sequence(organism).tsv :persist => true
                            protein_sequence.unnamed = true
                            protein_sequence.values_at(*self.ensembl)
                          end
  end

  property :sequence_length => :array2single do
    sequence.collect{|seq| seq.nil? ? nil : seq.length}
  end
end

