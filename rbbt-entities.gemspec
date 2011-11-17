# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{rbbt-entities}
  s.version = "1.0.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Miguel Vazquez"]
  s.date = %q{2011-11-17}
  s.description = %q{Entities for handling tsv files, caches, etc}
  s.email = %q{miguel.vazquez@cnio.es}
  s.extra_rdoc_files = [
    "LICENSE"
  ]
  s.files = [
    "LICENSE",
    "lib/rbbt/entity.rb",
    "lib/rbbt/entity/cnv.rb",
    "lib/rbbt/entity/gene.rb",
    "lib/rbbt/entity/genomic_mutation.rb",
    "lib/rbbt/entity/genotype.rb",
    "lib/rbbt/entity/misc.rb",
    "lib/rbbt/entity/mutated_isoform.rb",
    "lib/rbbt/entity/pmid.rb",
    "lib/rbbt/entity/protein.rb"
  ]
  s.homepage = %q{http://github.com/mikisvaz/rbbt-util}
  s.require_paths = ["lib"]
  s.rubygems_version = %q{1.6.2}
  s.summary = %q{Entities for the Ruby Bioinformatics Toolkit (rbbt)}
  s.test_files = ["test/test_helper.rb", "test/rbbt/entity/test_gene.rb", "test/rbbt/entity/test_genomic_mutation.rb", "test/rbbt/entity/test_mutated_isoform.rb", "test/rbbt/entity/test_protein.rb", "test/rbbt/test_entity.rb"]

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<rbbt-util>, [">= 0"])
      s.add_runtime_dependency(%q<rbbt-sources>, [">= 0"])
      s.add_runtime_dependency(%q<rbbt-dm>, [">= 0"])
    else
      s.add_dependency(%q<rbbt-util>, [">= 0"])
      s.add_dependency(%q<rbbt-sources>, [">= 0"])
      s.add_dependency(%q<rbbt-dm>, [">= 0"])
    end
  else
    s.add_dependency(%q<rbbt-util>, [">= 0"])
    s.add_dependency(%q<rbbt-sources>, [">= 0"])
    s.add_dependency(%q<rbbt-dm>, [">= 0"])
  end
end

