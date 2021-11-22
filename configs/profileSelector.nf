profiles {
	local {
		executor {
			name = 'local'
			cpus = 6
			//queuesize = 10
		}
	}

	chpc {
		executor {
			name = 'pbspro'
			queue = 'serial'
			queuesize = 10
		}

		singularity {
			enabled = true
			automounts = true
			cacheDir = "${HOME}/singularity"
		}

		params {
			containersDir = '/mnt/lustre/users/kesoh/containers/'
		}

		process {
			beforeScript = 'module load chpc/singularity/3.5.3'
			errorStrategy = { task.exitStatus in [135,255] ? 'retry' : 'finish' }
			maxErrors = '-1'
			maxRetries = 3
			clusterOptions = '-P CBBI1243 -l select=1 -m b'
			cpus = 24
			time = 90.m

			withLabel:ageEstimate {
				time = 90.m
				cpus = 24
				memory = 3.GB
			}

			withLabel:hapAnalysis {
				time = 10.m
				cpus = 24
				memory = 1.GB
			}

			withLabel:dmle {
				container = "${params.containersDir}dmle_latest.sif"
			}

                        withLabel:rcran {
                                container = "${params.containersDir}R.sif"
                        }
		}

	}
}

