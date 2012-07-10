require 'rbbt/entity'
require 'rbbt/entity/document'
require 'rbbt/sources/pubmed'

module PMID
  extend Entity
  include Document

  self.annotation :default_type

  self.format = "PMID"

  property :docid => :single do |*args|
    type = args.first || default_type
    ["PMID", self, type].compact * ":"
  end

  property :article => :array2single do
    PubMed.get_article(self).values_at(*self)
  end

  property :abstract => :array2single do
    article.collect{|a| a.nil? ? nil : a.abstract}
  end

  property :title => :array2single do
    article.collect{|a| a.nil? ? nil : a.title}
  end

  property :journal => :array2single do
    article.collect{|a| a.nil? ? nil : a.journal}
  end

  property :year => :array2single do
    article.collect{|a| a.nil? ? nil : a.year}
  end

  property :_get_text => :array2single do |*args|
    type = args.first || default_type

    case type.to_s
    when "full_text", 'fulltext'
      article.collect{|a| a.nil? ? nil : a.full_text}
    when "abstract"
      article.collect{|a| a.nil? ? nil : a.abstract }
    when "best"
      article.collect{|a| a.nil? ? nil : (a.full_text || a.text) }
    else
      article.collect{|a| a.nil? ? nil : a.text}
    end
  end

  property :pubmed_url => :single2array do
    "<a class='pmid' href='http://www.ncbi.nlm.nih.gov/pubmed/#{self}'>#{ self }</a>"
  end
  persist :pubmed_url 
end

