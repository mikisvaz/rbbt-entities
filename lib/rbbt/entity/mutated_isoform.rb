require 'rbbt/entity'
require 'rbbt/workflow'
require 'rbbt/sources/organism'
require 'rbbt/mutation/mutation_assessor'
require 'rbbt/entity/protein'
require 'rbbt/entity/gene'

module MutatedIsoform
  extend Entity
  self.annotation :organism

  self.format = "Mutated Isoform"

  property :protein do
    if Array === self
      Protein.setup(self.collect{|mutation| mutation.split(":").first}, "Ensembl Protein ID", organism)
    else
      Protein.setup(self.split(":").first, "Ensembl Protein ID", organism)
    end
  end

  property :change => :single2array do
    self.split(":").last
  end

  property :position => :single2array do
    if change.match(/[^\d](\d+)[^\d]/)
      $1.to_i
    else
      nil
    end
  end
  
  property :ensembl_protein_image_url => :single2array do
    ensembl_url = if organism == "Hsa" then "www.ensembl.org" else "#{organism.sub(/.*\//,'')}.archive.ensembl.org" end
    "http://#{ensembl_url}/Homo_sapiens/Component/Transcript/Web/TranslationImage?db=core;p=#{protein};_rmd=d2a8;export=svg"
  end

  ASTERISK = "*"[0]
  CONSECUENCES = %w(UTR SYNONYMOUS NOSTOP MISS-SENSE INDEL FRAMESHIFT NONSENSE)
  property :consecuence => :single2array do
    prot, change = self.split(":")

    case
    when change =~ /UTR/
      "UTR"
    when (change[0] == ASTERISK and not change[0] == change[-1])
      "NOSTOP"
    when (change[-1] == ASTERISK and not change[0] == change[-1])
      "NONSENSE"
    when change =~ /Indel/
      "INDEL"
    when change =~ /FrameShift/
      "FRAMESHIFT"
    when change[0] == change[-1]
      "SYNONYMOUS"
    else
      "MISS-SENSE"
    end
  end

  property :truncated => :array2single do
    @truncated ||= begin
                     protein2sequence_length = Misc.process_to_hash(self.protein.flatten){|list| list.sequence_length}
                     self.collect do |isoform_mutation|

                       next if isoform_mutation.consecuence != "FRAMESHIFT" and isoform_mutation.consecuence != "NONSENSE"
                       protein = isoform_mutation.protein
                       position = isoform_mutation.position
                       sequence_length = protein2sequence_length[protein]

                       case
                       when (sequence_length.nil? or position.nil?)
                         nil
                       when position < sequence_length.to_f * 0.7
                         true
                       else
                         false
                       end
                     end
                   end

  end

  property :mutation_assessor_scores => :array2single do
    @mutation_assesor_scores ||= begin
                                   missense = self.select{|mutation| mutation.consecuence == "MISS-SENSE"}

                                   correspondance = {}
                                   list = missense.zip(missense.protein.to "UniProt/SwissProt ID").collect do |mutation, uniprot|
                                     prot, change = mutation.split(":")
                                     next if uniprot.nil?
                                     uniprot_change = [uniprot, change]
                                     correspondance[uniprot_change] ||= []
                                     correspondance[uniprot_change] << mutation
                                     uniprot_change
                                   end.compact

                                   #return TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list) if list.empty?
                                   return [nil] * self.length if list.empty?

                                   tsv = MutationAssessor.chunked_predict(list.sort_by{|p| p * "_"})

                                   #return TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list) if tsv.nil? or tsv.empty? 
                                   return [nil] * self.length if tsv.empty?

                                   new = TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Func. Impact"], :type => :list)

                                   tsv.each do |key, values|
                                     correspondance[key.split(" ")].each do |mutation|
                                       new[mutation] = values["Func. Impact"]
                                     end
                                   end

                                   new.values_at *self
                                 end
  end

end


