#! /bin/env ruby


##################################################
# update: 20230220
# allows outputting only those within the operon

# update: 20230204
# multiple hits are considered


##################################################
require 'getoptlong'
require 'set'
require 'pry'
require 'colorize'

require 'util'
require 'Dir'
require 'SSW_bio'

$: << File.expand_path("/home-user/sswang/LHW-tools/extract/lib")
require 'util_extractCOG'


##################################################
def create_num_core_cutoff_range(str)
  errMsg("incorrect format for -n str. Should be in the format 2-10 or 2:10") if str !~ /[-:]/
  a1 = []
  a2 = str.split(/[-:]/)
  a1[0] = a2[0] =~ /^$/ ? 0 : a2[0].to_i
  a1[1] = a2[1] =~ /^$/ ? Float::INFINITY : a2[1].to_i
  range = Range.new(a1[0], a1[1])
  return(range)
end


def readCoreGeneFile(infile)
  h = Hash.new{|h,k|h[k]={}}
  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    gene = line_arr[0]
    if line_arr.size >= 3
      evalue, identity = line_arr[1,2].map{|i|i.to_f}
      h[gene][:evalue] = evalue
      h[gene][:identity] = identity
    else
      h[gene] = nil
    end
  end
  in_fh.close
  return(h)
end


def getCoor(locus, num_separated_genes=1)
  locus =~ /\d+$/
  coor = ($&.to_i)/num_separated_genes
  return(coor)
end


def findCoreGene(blast_file:, geneInfo:, c:, orgn:)
  subject2locus, locus2subject = Hash.new, Hash.new
  core_locus2subject = Hash.new
  core_subject2locus = Hash.new{|h,k|h[k]=[]}
  subject2locus = Hash.new{|h,k|h[k]=[]}
  #in_fh = File.open(blast_file, 'r')
  lines = File.readlines(blast_file)
  lines.each do |line|
    line.chomp!
    line_arr = line.split("\t")
    subject = line_arr[1]
    identity, evalue = line_arr.values_at(2,-2).map{|i|i.to_f}
    #if ! geneInfo[c].nil?
    if geneInfo.include?(c) and ! geneInfo[c].nil?
      next if evalue > geneInfo[c][:evalue]
      next if identity < geneInfo[c][:identity]
      core_subject2locus[c] << subject
      core_locus2subject[subject] = c
    end
    subject2locus[c] << subject
    locus2subject[subject] = c
  end
  #in_fh.close
  return([subject2locus, locus2subject, core_subject2locus, core_locus2subject])
end


def outputRes(core_sets:, subject2locus:, locus2subject:, dist_cutoff:, orgn:, orgn2seqfile:, seq_outdirs:, key_gene_2_abbr:, num_core_cutoff_range:, num_separated_genes:, include_genes:)
  # note there could be cases of >=2 sets of T3SS
  core_sets.each do |core_genes|
    next if core_genes.nil? or ! num_core_cutoff_range.include? core_genes.size
    #set = subject2locus.values.to_set.divide{|i,j| (getCoor(i,num_separated_genes)-getCoor(j,num_separated_genes)).abs <= dist_cutoff}.to_a.map{|i|i.to_a}
    #set = core_sets
    # genes = [[gene1, gene2], [gene3, gene4]]
    #geness = core_sets.sort_by{|i|i.to_a.size}.reverse
    geness = [core_genes]

    geness.each do |genes|
      if (genes & core_genes).size == core_genes.size
        if $is_print_orgn
          puts orgn
        elsif $is_print_gene
          #puts [orgn, include_genes.map{|gene| subject2locus.include?(gene) and ! (subject2locus[gene] & genes).empty? ? 1 : 0}].join("\t")
          puts [orgn, include_genes.select{|gene| subject2locus.include?(gene)}.map{|gene| !(subject2locus[gene] & genes).empty? ? 1 : 0}].join("\t")
        elsif $is_print_seq
          output_seq(genes:genes, orgn:orgn, locus2subject:locus2subject, orgn2seqfile:orgn2seqfile, seq_outdirs:seq_outdirs, key_gene_2_abbr:key_gene_2_abbr)
        else
          p [orgn, genes]
        end
      end
    end
  end
end


def get_orgn_2_seq_file(seq_files)
  h = Hash.new
  seq_files.map{|seq_file| h[getCorename(seq_file)] = seq_file}
  return(h)
end


def output_seq(genes:, orgn:, locus2subject:, orgn2seqfile:, seq_outdirs:, key_gene_2_abbr:)
  seqfile = orgn2seqfile[orgn]
  seq_objs = Hash.new
  seq_objs_ori = read_seq_file(seqfile)
  seq_objs_ori.each{|k,v| seq_objs[k.split(' ')[0]] = v } # there might be a space in seq title in the seq file like "brady|gene1 unknown"
  locus2subject.each do |locus, subject|
    next if not genes.include? locus if $is_only_within_operon
    if key_gene_2_abbr.include?(subject)
      outfile = File.join(seq_outdirs['key'], key_gene_2_abbr[subject]+'.fas')
      output_gene(outfile, seq_objs, orgn, locus)
    end
    outfile = File.join(seq_outdirs['all'], subject+'.fas')
    output_gene(outfile, seq_objs, orgn, locus)
  end
end


def output_gene(outfile, seq_objs, orgn, locus)
  out_fh = File.open(outfile, 'a')
  title = $is_print_locus ? locus : orgn
  out_fh.puts '>' + title + "\n" + seq_objs[locus].seq
  out_fh.close
end


def getAllQueries(geneInfo, type)
  #queries = `for i in \`find #{indir} -name '*blast8'\`; do cut -f1 $i; done \| sort \| uniq`.chomp
  genes = Array.new
  geneInfo.each_pair do |gene, v|
    genes << gene if gene =~ /^#{type}/
  end
  return(genes)
end


##################################################
indirs = Array.new
coreGeneFile = nil
type = nil
dist_cutoff = 20
seq_indirs = Array.new
seq_suffix = %w[protein gene fasta fas]
seq_files = Array.new
species_included = Array.new
key_gene_2_abbr = Hash.new
num_core_cutoff_range = Range.new(2,Float::INFINITY)

seq_outdir = nil
seq_outdirs = Hash.new
is_force = false

$is_print_orgn = false
$is_print_gene = false
$is_print_seq = false
$is_print_locus = false
$is_only_within_operon = false

geneInfo = Hash.new
orgn2seqfile = Hash.new


##################################################
argv = Marshal.load(Marshal.dump(ARGV))

opts = GetoptLong.new(
  ['--indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--core_gene', GetoptLong::REQUIRED_ARGUMENT],
  ['-t', '--type', GetoptLong::REQUIRED_ARGUMENT],
  ['-d', '--dist', GetoptLong::REQUIRED_ARGUMENT],
  ['-n', '--num_core_gene', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--seq_outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--include_list', '--species_include_list', GetoptLong::REQUIRED_ARGUMENT],
  ['--key_gene_list', GetoptLong::REQUIRED_ARGUMENT],
  ['--print_orgn', GetoptLong::NO_ARGUMENT],
  ['--print_gene', GetoptLong::NO_ARGUMENT],
  ['--print_seq', GetoptLong::NO_ARGUMENT],
  ['--print_locus', GetoptLong::NO_ARGUMENT],
  ['--only_within_operon', '--within_operon', '--only_operon', '--operon', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--indir'
      indirs << value.split(',')
      indirs.flatten!
    when '--core_gene'
      coreGeneFile = value
    when '-t', '--type'
      type = value
    when '-d', '--dist'
      dist_cutoff = value.to_i
    when '-n', '--num_core_gene'
      num_core_cutoff_range = create_num_core_cutoff_range(value)
    when '--seq_indir'
      seq_indirs << value.split(',')
      seq_indirs.flatten!
    when '--seq_outdir'
      seq_outdir = value
    when '--force'
      is_force = true
    when '--include_list', '--species_include_list'
      species_included = read_list(value).keys
    when '--key_gene_list'
      key_gene_2_abbr = readTbl(value)[0]
    when '--print_orgn'
      $is_print_orgn = true
    when '--print_gene'
      $is_print_gene = true
    when '--print_seq'
      $is_print_seq = true
    when '--print_locus'
      $is_print_locus = true
    when '--only_within_operon', '--within_operon', '--only_operon', '--operon'
      $is_only_within_operon = true
  end
end


##################################################
if $is_print_seq
  if seq_indirs.empty? or seq_outdir.nil?
    raise "--seq_indir and --seq_outdir have to be provided if --print_seq is specified."
  end
end


##################################################
geneInfo = readCoreGeneFile(coreGeneFile)
#geneInfo = gene2evalue.merge(gene2evalue){|k, v| v = v.to_f if not v.nil?}

if not seq_indirs.empty? and not seq_outdir.nil?
  seq_indirs.each do |seq_indir|
    seq_files << read_infiles(seq_indir, seq_suffix)
    seq_files.flatten!
  end
  orgn2seqfile = get_orgn_2_seq_file(seq_files)
  mkdir_with_force(seq_outdir, is_force) # create seq_outdir

  argv_outfile = File.join(seq_outdir, 'argv')
  argv_outfh = File.open(argv_outfile, 'w')
  argv_outfh.puts argv.join(' ')
  argv_outfh.close

  seq_outdirs['all'] = File.join(seq_outdir, 'all')
  seq_outdirs['key'] = File.join(seq_outdir, 'key')
  seq_outdirs.map{|k, v| mkdir_with_force(v, is_force) }
end

unless species_included.empty?
  $out2in, $in2out = getSpeciesNameRelaWithoutFile(species_included)
else
  $out2in, $in2out = getSpeciesNameRelaWithoutFile(seq_files.map{|i|[getCorename(i), 1]}.to_h)
end


##################################################
if $is_print_gene
  include_genes = getAllQueries(geneInfo, type)
  puts ['', include_genes].flatten.join("\t")
end


##################################################
indirs.each do |indir|
  Dir.foreach(indir) do |orgn|
    next if orgn =~ /^\./
    sub_indir = File.join(indir, orgn)
    subject2locus, locus2subject, core_subject2locus, core_locus2subject = Hash.new, Hash.new, Hash.new, Hash.new
    Dir.foreach(sub_indir) do |b|
      next if b =~ /^\./
      c = getCorename(b)
      blast_file = File.join(sub_indir, b)
      h1, h2, h3, h4 = findCoreGene(blast_file:blast_file, geneInfo:geneInfo, c:c, orgn:orgn)
      subject2locus.merge!(h1)
      locus2subject.merge!(h2)
      core_subject2locus.merge!(h3)
      core_locus2subject.merge!(h4)
    end

    # "num_separated_genes" is the num of genes separated, e.g., some genome annotation use five as the unit in locus name btwn genes
    num_separated_genes = locus2subject.keys.all?{|locus|locus =~ /[05]$/} ? 5 : 1


    core_sets = core_subject2locus.values.flatten.to_set.divide{|i,j| (getCoor(i,num_separated_genes)-getCoor(j,num_separated_genes)).abs <= dist_cutoff}.to_a.map{|i|i.to_a}
    core_sets = core_sets.sort_by{|i|i.size}.reverse
    #core_sets.flatten!(1)
    
    outputRes(core_sets:core_sets, subject2locus:subject2locus, locus2subject:locus2subject, dist_cutoff:dist_cutoff, orgn:orgn, orgn2seqfile:orgn2seqfile, seq_outdirs:seq_outdirs, key_gene_2_abbr:key_gene_2_abbr, num_core_cutoff_range:num_core_cutoff_range, num_separated_genes:num_separated_genes, include_genes:include_genes)
  end
end


