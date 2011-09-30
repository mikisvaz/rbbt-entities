require 'rbbt/annotations'

module Entity
  class << self
    attr_accessor :formats
  end
  self.formats = {}
  
  def self.extended(base)
    base.extend Annotation unless Annotation === base

    Entity.formats[base.to_s] = base
    base.module_eval do
      class << self
        attr_accessor :template
        alias prev_entity_extended extended
      end

      def self.extended(data)
        prev_entity_extended(data)
        data.extend AnnotatedArray
      end

      def self.format=(formats)
        formats = [formats] unless Array === formats
        formats.each do |format|
          Entity.formats[format] = self
        end
      end
    end
  end
end

