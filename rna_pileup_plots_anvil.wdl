version 1.0
workflow samtools_stats {
    input{
        File bam
        Int memory
        String chrom
        String tissue
        Int start
        Int end
        Int num_threads
        Int disk_space
        Int num_preempt
    }
	
    call samtools_stats{
        input:
            bam = bam,
            memory = memory,
            chrom = chrom,
            tissue = tissue,
            start = start,
            end = end,
            num_threads = num_threads,
            disk_space = disk_space,
            num_preempt = num_preempt
    }
    output {
    	File stats = samtools_stats.stats
    }
}
task samtools_stats {
    input{
    File bam
    Int memory
    String chrom
    String tissue
    Int start
    Int end
    Int num_threads
    Int disk_space
    Int num_preempt
    }    
    command<<<
        samtools index ~{bam}
        samtools depth -a -r ~{chrom}:~{start}-~{end} ~{bam} > ~{tissue}.txt
    >>>

    output {
        File stats = "${tissue}.txt"
    }

    runtime {
        docker: "mschatz/wga-essentials"
        cpu: "${num_threads}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Bohan Ni"
    }
}