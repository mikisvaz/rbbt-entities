#!/usr/bin/env ruby
require 'rbbt/util/open'

file = ARGV.shift

Open.read(file).each do |line|
  next if line =~ /^#/

  chr, pos, id, ref, mut, score = line.split(/\t/)

  mut = mut + '-' * (ref.length - mut.length) if ref.length > mut.length

  puts [chr, pos, mut, score] * ":"
end
