require 'rbbt/entity'
require 'rbbt/resource/path'
require 'rbbt/entity/genotype'
require 'rbbt/workflow/rest/render'
require 'rbbt/sources/barcode'

module Study
  extend Entity
  include RbbtHTMLHelpers

  self.annotation :workflow

  class << self
    alias prev_study_extended extended
  end

  def self.load(study, workflow = nil)
    study_dir = case
                when (not defined?(Rbbt))
                  File.join(ENV["HOME"], '.studies')
                when Rbbt.etc.study_dir.exists?
                  Rbbt.etc.study_dir.read.chomp
                else
                  Rbbt.studies.find
                end
    
    path = Study.setup(File.join(study_dir, study), workflow)
    path.extend Path
    raise "Study '#{ study }' not found in '#{ path }'" unless path.exists?

    path
  end

  def self.extended(object)
    self.send(:prev_study_extended, object)
    setup_file = File.join(object, 'setup.rb')
    if File.exists? setup_file
      class << object; self; end.class_eval Open.read(setup_file)
    end
    object
  end

  def name
    File.basename self
  end

  def metadata_file
    File.join(self, 'metadata.yaml')
  end

  def genotype_files
    Dir.glob(File.join(self, 'genotypes', '*')).reject{|p| p[0].chr == '.'}
  end

  def metadata
    if File.exists? metadata_file
      File.open(metadata_file){|f|
        YAML.load(f)
      }
    else
      {}
    end
  end

  def expressed_genes(threshold = nil)
    threshold ||= 0.85
    @expressed_genes ||= {}
    @expressed_genes[threshold] ||= if metadata.include?(:barcode_tissue) and not metadata[:barcode_tissue].nil?
                                      tissues = metadata[:barcode_tissue]
                                      tissues = [tissues] unless Array === tissues
                                      transcriptome = Barcode.transcriptome.tsv(:persist => true, :type => :list, :cast => :to_f, :unnamed => true)
                                      tissue_pos = tissues.collect{|tissue| transcriptome.identify_field tissue}
                                      transcriptome.select{|probe, values| values.values_at(*tissue_pos).select{|value| value >= threshold}.any? }.keys.ensembl
                                    else
                                      nil
                                    end
  end

  def cohort
    organism = metadata[:organism]
    watson = metadata[:watson]
    @cohort ||= genotype_files.collect do |f| 
      name = File.basename(f)
      genomic_mutations = Open.read(f).split("\n").reject{|l| l.empty? }
      GenomicMutation.setup(genomic_mutations, name, organism, watson)
    end.tap{|cohort| cohort.extend Genotype::Cohort}
  end

  def affected_genes
    cohort.metagenotype.genes.flatten.compact.uniq
  end

  def damaged_genes(*args)
    cohort.damaged_genes(*args)
  end

  def recurrent_genes
    cohort.recurrent_genes.annotate(cohort.recurrent_genes + Misc.counts(cohort.metagenotype.genes.compact.flatten).select{|g,c| c > 1}.collect{|g,c| g})
  end

  def relevant_genes(options = {})
    damaged_genes = self.damaged_genes(options)
    recurrent_genes = self.recurrent_genes
    cancer_related_genes = self.affected_genes.reject{|g| g.related_cancers.nil? or g.related_cancers.empty?}
    expressed_genes = self.expressed_genes(options[:expression_threshold])


    if expressed_genes.nil?
      damaged_genes.annotate (damaged_genes + recurrent_genes + cancer_related_genes).uniq
    else
      damaged_genes.annotate(damaged_genes.subset(expressed_genes) + recurrent_genes + cancer_related_genes).uniq
    end
  end

  def table(list, *fields)
    fields = fields.first if fields.length == 0 and Array === fields.first 
    rows = list.subset(affected_genes).collect{|entity|
      [entity.id, entity.tsv_values(fields)]
    }

    workflow_partial('partials/_table', workflow, :rows => rows, :header => fields)
  end

end
