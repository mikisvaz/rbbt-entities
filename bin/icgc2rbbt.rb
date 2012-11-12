#!/usr/bin/env ruby
require 'rbbt/util/open'

file = ARGV.shift
directory = ARGV.shift

genotypes = {}
Open.read(file).split("\n").each do |line|
  next if line =~ /^Cancer Type/

  chr, pos, ref, mut, sample = line.split(/\t/).values_at 2, 3, 6, 10, 28 

  chr.sub!(/chr/,'')
  mut = '-' * (mut.length - 1) if mut =~/^-[ACGT]/

  genotypes[sample] ||= []
  genotypes[sample] << [chr, pos, mut] * ":"
end

genotypes.each do |sample, mutations|
  mutations.uniq!
  Open.write(File.join(directory, sample), mutations.uniq * "\n")
end
