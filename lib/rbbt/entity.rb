require 'rbbt/annotations'

module Entity
  class << self
    attr_accessor :formats, :entity_property_cache, :entity_list_cache
  end

  self.entity_property_cache = "cache/entity_property"
  self.entity_list_cache = "cache/entity_list"
  self.formats = {}
  
  UNPERSISTED_PREFIX = "entity_unpersisted_property_"

  def self.extended(base)
    base.extend Annotation unless Annotation === base

    Entity.formats[base.to_s] = base
    base.module_eval do
      if not methods.include? "prev_entity_extended"
        class << self
          attr_accessor :template, :list_template, :action_template, :list_action_template, :keep_id
          alias prev_entity_extended extended 
        end 

        def self.extended(data)
          prev_entity_extended(data)

          if Array === data and 
            not AnnotatedArray === data and 
            not (Annotated === data.compact.first and (data.annotation_types - data.compact.first.annotation_types).any?)

            data.extend AnnotatedArray
          end

          data
        end
      end

      def self.format=(formats)
        formats = [formats] unless Array === formats
        formats.each do |format|
          Entity.formats[format] = self
        end
      end

      def clean_annotations
        case
        when self.nil?
          self
        when String === self
          "" << self
        when self.respond_to?(:annotated_array_clean_collect)
          self.annotated_array_clean_collect{|e| e.respond_to?(:clean_annotations)? e.clean_annotations : e}
        when Array === self
          self.collect{|e| e.respond_to?(:clean_annotations)? e.clean_annotations : e}
        when (Fixnum === self)
          0 + self
        when (Float === self)
          0.0 + self
        else
          raise "Unknown casuistic in clean_annotations for object: #{self.inspect}"
        end
      end

      def to_yaml(*args)
        clean_annotations.to_yaml(*args)
      end

      def marshal_dump
        clean_annotations
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

      def self.property(name, &block)
        case
        when (Hash === name and name.size == 1)
          name, type = name.collect.first
        when (String === name or Symbol === name)
          type = :single
        else
          raise "Format of name ( => type) not understood: #{name.inspect}"
        end

        name = name.to_s unless String === name

        persisted_name = UNPERSISTED_PREFIX + name.to_s
        self.remove_method persisted_name if methods.include? persisted_name

        case type
        when :single, :single2array
          single_name = "_single_" << name
          define_method single_name, &block 
          define_method name do |*args|
            if Array === self
              self.collect{|e| e.send(name, *args)}
            else
              self.send(single_name, *args)
            end
          end
        when :array, :array2single
          ary_name = "_ary_" << name
          define_method ary_name, &block 

          define_method name do |*args|
            case
            when Array === self
              self.send(ary_name, *args)
            when (Array === self.container and self.container.respond_to? ary_name)
              res = self.container.send(name, *args)
              if Hash === res
                res[self]
              else
                res[self.container_index]
              end
            else
              res = self.make_list.send(ary_name, *args)
              Hash === res ? res[self] : res[0]
            end
          end

        end
      end

      def self.persist(method_name, type = nil, options = {})
        type = :memory if type.nil?
        options = Misc.add_defaults options, :dir => Entity.entity_property_cache

        orig_name = UNPERSISTED_PREFIX + method_name.to_s
        alias_method orig_name, method_name unless instance_methods.include? orig_name

        define_method method_name do |*args|
          persist_name = __method__.to_s << ":" << self.id

          persist_options = options
          persist_options = persist_options.merge(:other => {:args => args}) if args.any?

          Persist.persist(persist_name, type, persist_options) do
            self.send(orig_name, *args)
          end
        end
      end

    end 
  end
end

