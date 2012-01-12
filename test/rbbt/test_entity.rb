require File.expand_path(File.dirname(__FILE__) + '/../test_helper')
require 'rbbt/entity'
require 'rbbt/util/tmpfile'
require 'test/unit'

module ReversableString
  extend Entity
  
  property :reverse_text_ary => :array do
    $count += 1
    self.collect{|s| s.reverse}
  end
 
  property :reverse_text_single => :single do
    $count += 1
    self.reverse
  end

  property :reverse_text_ary_p => :array do
    $count += 1
    self.collect{|s| s.reverse}
  end
 
  property :reverse_text_single_p => :single do
    $count += 1
    self.reverse
  end

  property :reverse_text_ary_p_array => :array do
    $count += 1
    self.collect{|s| s.reverse}
  end
 
  persist :reverse_text_ary_p
  persist :reverse_text_ary_p
  persist :reverse_text_ary_p
  persist :reverse_text_single_p

  persist :reverse_text_ary_p_array, :array, :dir => TmpFile.tmp_file
end
class TestEntity < Test::Unit::TestCase

  def test_property_ary
    a = ["String1", "String2"]
    a.extend ReversableString

    $count = 0

    assert_equal "2gnirtS", a.reverse_text_ary.last
    assert_equal 1, $count
    assert_equal "2gnirtS", a[1].reverse_text_ary
    assert_equal 2, $count
  end

  def test_property_single
    a = ["String1", "String2"]
    a.extend ReversableString

    $count = 0

    assert_equal "2gnirtS", a.reverse_text_single.last
    assert_equal 2, $count
    assert_equal "2gnirtS", a[1].reverse_text_single
    assert_equal 3, $count
  end

  def test_property_ary_p
    a = ["String1", "String2"]
    a.extend ReversableString

    $count = 0

    assert_equal "2gnirtS", a.reverse_text_ary_p.last
    assert_equal "2gnirtS", a[1].reverse_text_ary_p
    assert_equal 1, $count
  end

  def test_property_single_p
    a = ["String1", "String2"]
    a.extend ReversableString

    $count = 0

    assert_equal "2gnirtS", a.reverse_text_single_p.last
    assert_equal 2, $count
    assert_equal "2gnirtS", a[1].reverse_text_single_p
    assert_equal 2, $count
  end

  def test_property_ary_p_array
    a = ["String1", "String2"]
    a.extend ReversableString

    $count = 0

    assert_equal "2gnirtS", a.reverse_text_ary_p_array.last
    assert_equal 1, $count
    assert_equal "2gnirtS", a.reverse_text_ary_p_array.last
    assert_equal 1, $count
  end
end
