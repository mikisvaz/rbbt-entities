require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/genomic_mutation'

class TestGenomicMutation < Test::Unit::TestCase
  MUTATION = GenomicMutation.setup("10:124745844:A:158", "Test", "Hsa/jun2011")
  SPLICING = GenomicMutation.setup("18:14787040:A", "Test", "Hsa/jun2011")
  GENOTYPE = GenomicMutation.setup(Rbbt.data.genotype.list, "Test", "Hsa/jun2011")

  def test_genes
    assert GENOTYPE.genes.flatten.to("Associated Gene Name").include? "PSTK"
  end

  def test_consolidate
    assert GENOTYPE.genes.consolidate.to("Associated Gene Name").include? "PSTK"
  end

  def test_mutated_isoforms
    assert MUTATION.mutated_isoforms.length > 1
    assert ["PSTK"], MUTATION.mutated_isoforms.protein.gene.to("Associated Gene Name").uniq
  end

  def test_exon_junction
    assert(!(MUTATION.in_exon_junction?))
    assert SPLICING.in_exon_junction?
  end

  def test_over_gene
    assert MUTATION.over_gene? Gene.setup("PSTK", "Associated Gene Name", "Hsa/jun2011").ensembl
    assert(!(SPLICING.over_gene? Gene.setup("PSTK", "Associated Gene Name", "Hsa/jun2011").ensembl))
  end

end


