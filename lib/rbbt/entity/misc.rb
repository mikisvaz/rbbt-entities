require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/go'
require 'rbbt/sources/organism'
require 'rbbt/entity/gene'

module GOTerm
  extend Entity
  self.annotation :organism

  self.format = ["GO Term", "GO ID"]

  def name
    if Array === self
      self.collect{|id| GO.id2name(id)}
    else
      GO.id2name(self)
    end
  end

  def genes
    go2genes = Organism.gene_go(organism).tsv(:key_field => "GO ID", :fields => ["Ensembl Gene ID"], :merge => true, :persist => true)
    go2genes.unnamed = true
    Gene.setup(go2genes[self].first, "Ensembl Gene ID", organism)
  end
end
