workflow genomicsVarTrix {
  ## required files
  File cell_barcodes
  File bamFile
  File baiFile

  File vcfFile
  File vcfiFile

  ## reference files
  File refFasta
  File refFastai


  String sample_id

  ## runtime
  String docker
  Int preemptible_tries
  Int disk_size

call VarTrix {
    input: 
      cell_barcodes=cell_barcodes,
      bamFile=bamFile,
      baiFile=baiFile,
      RefFasta=refFasta, 
      RefIndex=refFastai,
      vcfFile=vcfFile,
      vcfiFile=vcfiFile,
      sample_name=sample_id,
      docker=docker,
      preemptible_tries=preemptible_tries,
      disk_size=disk_size
  }
}

task VarTrix {

  File cell_barcodes
  File bamFile
  File baiFile
  File RefFasta
  File RefIndex
  File vcfFile
  File vcfiFile

  String sample_name
  String docker

  Int preemptible_tries
  Int disk_size

  command <<<
    vartrix -v ${vcfFile} 
              -b ${bamFile} 
              -f ${RefFasta} 
              -c ${cell_barcodes}
              --umi
              -s alt_frac
              -o ${sample_name}

    vawk '{print $1,$2}' vcfFile > SNV.loci.txt
    sed -i 's/\s/:/g' SNV.loci.txt 

  >>>

  output {
    File output_file = "${sample_name}"
    File snv_loci = "SNV.loci.txt"
  }
    runtime {
        preemptible: preemptible_tries
        memory: "20 GB"
        cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
      }
}




