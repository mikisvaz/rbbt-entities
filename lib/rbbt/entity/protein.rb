require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/statistics/hypergeometric'
require 'rbbt/network/paths'

Workflow.require_workflow "Translation"

module Protein
  extend Entity
  include Entity::Adjacent
  include Entity::Enriched

  self.annotation :format
  self.annotation :organism

  def gene
    Gene.setup(to("Ensembl Protein ID"), "Ensembl Protein ID", organism)
  end

  def to(new_format)
    return self if format == new_format
    if Array === self
      Protein.setup(Translation.job(:translate_protein, "", :organism => organism, :proteins => self, :format => new_format).exec, new_format, organism)
    else
      Protein.setup(Translation.job(:translate_protein, "", :organism => organism, :proteins => [self], :format => new_format).exec.first, new_format, organism)
    end
  end
end

