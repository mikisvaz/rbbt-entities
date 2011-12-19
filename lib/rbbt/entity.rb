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
        data.extend AnnotatedArray if Array === data
        data
      end

      def self.format=(formats)
        formats = [formats] unless Array === formats
        formats.each do |format|
          Entity.formats[format] = self
        end
      end

      def clean_annotations
        case
        when Array === self
          self.annotated_array_clean_collect{|e| e.respond_to?(:clean_annotations)? e.clean_annotations : e}
        when String === self
          "" << self
        end
      end

      def consolidate
        self.inject(nil){|acc,e| 
          if acc.nil?
            acc = e
          else
            acc.concat e
          end
        }
      end
    end
  end

  def property(name, &block)
    case
    when (Hash === name and name.size == 1)
      name, type = name.collect.first
    when (String === name or Symbol === name)
      type = :both
    else
      raise "Format of name ( => type) not understood: #{name.inspect}"
    end

    name = name.to_s unless String === name

    case type
    when :both
      self.module_eval do define_method name, &block end
    when :array
      self.module_eval do 
        ary_name = "_ary_" << name
        define_method ary_name, &block 
        define_method name do |*args|
          raise "Method #{ name } only defined for array" unless Array === self
          self.send(ary_name, *args)
        end
      end
    when :single
      self.module_eval do 
        single_name = "_single_" << name
        define_method single_name, &block 
        define_method name do |*args|
          raise "Method #{ name } not defined for array" if Array === self
          self.send(single_name, *args)
        end
      end
    when :single2array
      self.module_eval do 
        single_name = "_single_" << name
        define_method single_name, &block 
        define_method name do |*args|
          if Array === self
            collect{|e| e.send(single_name, *args)}
          else
            self.send(single_name, *args)
          end
        end
      end
    when :array2single
      self.module_eval do 
        ary_name = "_ary_" << name
        define_method ary_name, &block 
        define_method name do |*args|
          case
          when Array === self
            self.send(ary_name, *args)
          when (Array === self.container and self.container.respond_to? ary_name)
            res = self.container.send(ary_name, *args)
            if Hash === res
              res[self]
            else
              pos = self.container.index self
              res[pos]
            end
          else
            res = self.make_list.send(ary_name, *args)
            Hash === res ? res[self] : res[0]
          end
        end
      end
    end
  end
end

