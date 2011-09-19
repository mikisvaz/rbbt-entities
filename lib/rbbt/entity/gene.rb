require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/statistics/hypergeometric'
require 'rbbt/network/paths'

Workflow.require_workflow "Translation"

module Gene
  extend Entity
  include Entity::Adjacent
  include Entity::Enriched

  self.annotation :format
  self.annotation :organism

  self.format = Organism::Hsa.identifiers.all_fields

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
      to!(new_format).collect{|v| v.nil? ? nil : v.first}
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
end
