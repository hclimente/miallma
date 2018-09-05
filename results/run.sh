module add Java

../scripts/0.1_prepare_expression.nf -resume -profile quickie
../scripts/0.2_get_drivers.nf -resume -profile quickie
../scripts/1_run_spada.nf -resume -profile cluster
