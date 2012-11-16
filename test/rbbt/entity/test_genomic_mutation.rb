require File.expand_path(File.dirname(__FILE__) + '/../../test_helper')

require 'test/unit'
require 'rbbt/util/tmpfile'
require 'test/unit'
require 'rbbt/entity/genomic_mutation'

class TestGenomicMutation < Test::Unit::TestCase
  MUTATION = GenomicMutation.setup("10:124745844:A:158", "Test", "Hsa/jun2011")
  NON_CODING_MUTATION = GenomicMutation.setup("12:52844451:T", "Test", "Hsa/jun2011")
  SPLICING = GenomicMutation.setup("18:14787040:A", "Test", "Hsa/jun2011")

  FORWARD_STRAND_FIRST_EXON_START = GenomicMutation.setup("10:89622870", "Test", "Hsa/jun2011")
  FORWARD_STRAND_FIRST_EXON_END   = GenomicMutation.setup("10:89624305", "Test", "Hsa/jun2011")
  FORWARD_STRAND_LAST_EXON_START  = GenomicMutation.setup("10:89725044", "Test", "Hsa/jun2011")
  FORWARD_STRAND_LAST_EXON_END    = GenomicMutation.setup("10:89731687", "Test", "Hsa/jun2011")

  REVERSE_STRAND_FIRST_EXON_START = GenomicMutation.setup("2:198256698", "Test", "Hsa/jun2011")
  REVERSE_STRAND_FIRST_EXON_END   = GenomicMutation.setup("2:198257185", "Test", "Hsa/jun2011")
  REVERSE_STRAND_LAST_EXON_START  = GenomicMutation.setup("2:198299696", "Test", "Hsa/jun2011")
  REVERSE_STRAND_LAST_EXON_END    = GenomicMutation.setup("2:198299815", "Test", "Hsa/jun2011")


  def test_mutated_isoforms
    assert MUTATION.mutated_isoforms.length > 1
    assert_equal ["PSTK"], MUTATION.mutated_isoforms.protein.gene.to("Associated Gene Name").uniq
  end

  def test_exon_junction
    assert(!(MUTATION.in_exon_junction?))
    assert SPLICING.in_exon_junction?

    assert(!(FORWARD_STRAND_FIRST_EXON_START.in_exon_junction?))
    assert( (FORWARD_STRAND_FIRST_EXON_END.in_exon_junction?))
    assert( (FORWARD_STRAND_LAST_EXON_START.in_exon_junction?))
    assert(!(FORWARD_STRAND_LAST_EXON_END.in_exon_junction?))

    assert(!(REVERSE_STRAND_FIRST_EXON_START.in_exon_junction?))
    assert( (REVERSE_STRAND_FIRST_EXON_END.in_exon_junction?))
    assert( (REVERSE_STRAND_LAST_EXON_START.in_exon_junction?))
    assert(!(REVERSE_STRAND_LAST_EXON_END.in_exon_junction?))
  end

  def test_over_gene
    assert MUTATION.over_gene? Gene.setup("PSTK", "Associated Gene Name", "Hsa/jun2011").ensembl
    assert(!(SPLICING.over_gene? Gene.setup("PSTK", "Associated Gene Name", "Hsa/jun2011").ensembl))
  end

  def test_type
    reference = "A"
    mutation = GenomicMutation.setup("1:1727802:A", "Test", "Hsa/may2009", true)
    assert_equal 'none', mutation.type
    mutation = GenomicMutation.setup("1:1727802:C", "Test", "Hsa/may2009", true)
    assert_equal 'transversion', mutation.type
    mutation = GenomicMutation.setup("1:1727802:T", "Test", "Hsa/may2009", true)
    assert_equal 'transversion', mutation.type
    mutation = GenomicMutation.setup("1:1727802:G", "Test", "Hsa/may2009", true)
    assert_equal 'transition', mutation.type
  end

  def test_coding
    assert MUTATION.coding?
    ddd NON_CODING_MUTATION.affected_exons
    assert(!NON_CODING_MUTATION.coding?)
  end

end


