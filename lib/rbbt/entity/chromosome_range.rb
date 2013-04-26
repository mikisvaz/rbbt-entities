require 'rbbt'
require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/entity/gene'

Workflow.require_workflow "Sequence"

module ChromosomeRange
  extend Entity

  self.annotation :organism
 
  self.format = "Chromosome Range"


  def self.text_to_unit(text)
    text = text.sub('^', '+')
    base = text.to_f
    case
    when text =~ /KB?$/i
      base * 1000
    when text =~ /MB?$/i
      base * 1000_000
    when text =~ /^\d+(\.\d+)?(e\+\d+)?$/
      base
    else 
      raise "Text format not understood: #{ text }"
    end.to_i
  end

  property :unit => :array2single do
    res = self.collect{|range|
      chr, start, eend, *rest = range.split(":")
      ([chr, ChromosomeRange.text_to_unit(start), ChromosomeRange.text_to_unit(eend)].concat rest) * ":"
    }
    self.annotate res
    res.extend AnnotatedArray if AnnotatedArray === self
    res
  end

  property :genes => :array2single do
    Sequence.job(:genes_at_genomic_ranges, "ChromosomeRange", :organism => organism, :ranges => self.unit).run.tap{|t| t.namespace = organism}.chunked_values_at self.unit
  end

  property :ensembl_browser => :single2array do
    "http://#{Misc.ensembl_server(self.organism)}/#{Organism.scientific_name(organism).sub(" ", "_")}/Location/View?db=core&r=#{chromosome}:#{start}-#{eend}"
  end

  property :ucsc_browser => :single2array do
    "http://genome.ucsc.edu/cgi-bin/hgTracks?db=#{Organism.hg_build(organism)}&position=chr#{chromosome}:#{start}-#{eend}"
  end

  property :chromosome => :array2single do
    self.clean_annotations.collect{|r| r.split(":")[0]}
  end

  property :start => :array2single do
    self.clean_annotations.collect{|r| r.split(":")[1].to_i}
  end

  property :begin => :array2single do
    start
  end

  property :eend => :array2single do
    self.clean_annotations.collect{|r| r.split(":")[2].to_i}
  end

  property :end => :array2single do
    self.eend
  end

  property :range => :array2single do
    start.zip(self.end).collect{|s,e| (s..e)}
  end

  property :over_chromosome_range? => :array2single do |r|
    range_chr, start, eend = r.split(":")
    range_range = (start.to_i..eend.to_i)
    chromosome.zip(range).collect do |chr,range| 
      chr == range_chr and ( range.include?(range_range.begin) or range.include?(range_range.end) or range_range.include?(range.begin) or range_range.include?(range.end))
    end
  end

end
