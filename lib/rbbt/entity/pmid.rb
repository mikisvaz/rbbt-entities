require 'rbbt/entity'
require 'rbbt/entity/document'
require 'rbbt/sources/pubmed'
require 'rbbt/sources/gscholar'

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
    PubMed.get_article(self).chunked_values_at(self)
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

  property :cites => :single2array do
    if title
      begin
        GoogleScholar.number_cites(title)
      rescue
        nil
      end
    else
      nil
    end
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

  property :bibtex => :array2single do
    PubMed.get_article(self).chunked_values_at(self).collect do |article|
      article.bibtex
    end
  end

end
