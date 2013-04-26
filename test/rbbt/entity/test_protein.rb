require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'rbbt'
require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/protein'
require 'rbbt/entity/mutated_isoform'

class TestProtein < Test::Unit::TestCase
  PROTEIN = Protein.setup("ENSP00000275493", "Ensembl Protein ID", "Hsa/jun2011")
  PROTEIN_ARRAY = Protein.setup(["ENSP00000275493"], "Ensembl Protein ID", "Hsa/jun2011")

  def test_clean_annotations
    assert Protein === PROTEIN
    assert(!(Protein === PROTEIN.clean_annotations))
    assert Gene === PROTEIN.gene
    assert(!(Protein === PROTEIN.gene))
  end

  def test_gene
    assert_equal "EGFR", PROTEIN.gene.to("Associated Gene Name")
  end

  def test_protein
    assert MutatedIsoform.setup(["ENSP00000322422:K168FrameShift"], "Hsa/may2009").info.include?(:organism)
    assert_equal ["ENSP00000322422"], MutatedIsoform.setup(["ENSP00000322422:K168FrameShift"], "Hsa/may2009").protein
  end

  def test_ortholog
    assert_equal "ENSMUSP00000020329", PROTEIN.ortholog("Mmu/jun2011")
  end


end


