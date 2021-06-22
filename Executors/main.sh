
    #!/bin/sh

    #obtain id of preprocess job
    cd /home/emiliano/projects/def-cdesouza/Lab/MRNA_pipe/Pipe/rnaSeq/code
    git checkout ${BRANCH}
    git pull --all 

    if [ "$PREPROCESS" == "true" ]
    then 
        ###################### PREPROCESS ######################
        Preprocess__id=$(sbatch --parsable --job-name=preprocess_${DATA_NAME} \
        --output=../../rnaSeq/output/preprocess/${DATA_NAME}/preprocess-${DATA_NAME}.out \
        --export=DATASET=${DATA_NAME} ../../rnaSeq/code/preprocess/preprocess_job.sh) 

        echo "$Preprocess__id"
        start_time=$(date +%M.%S)
        while [ $(sacct -j ${Preprocess__id} --format=State| tail -n1 | xargs) != "COMPLETED" ]
        do 

            sleep 30s
            end_time=$(date +%M.%S)
            elapsed=$(echo "scale=0; $end_time - $start_time" | bc)

            echo "Job Pending time $elapsed "
        done 

        echo "====================== PREPROCESS CONSOLE LOG ======================  "

        cat "../../rnaSeq/output/preprocess/${DATA_NAME}/preprocess-${DATA_NAME}.out"
    fi

    if [ "$CLUSTER" == "true" ]
    then
        ###################### CLUSTER ######################
        echo "====================== CLUSTER ANALYSIS ======================  "

        cluster1_id=$(sbatch --parsable --job-name=sc3_${DATA_NAME} \
        --output=../../rnaSeq/output/cluster/sc3/${DATA_NAME}/sc3-${DATA_NAME}.out \
        --export=DATASET=${DATA_NAME},CLUSTERMETHOD=sc3 ../../rnaSeq/code/cluster/cluster_job.sh) 

        echo "sc3 job id $cluster1_id"

        cluster2_id=$(sbatch --parsable --job-name=seurat_${DATA_NAME} \
        --output=../../rnaSeq/output/cluster/seurat/${DATA_NAME}/seurat-${DATA_NAME}.out \
        --export=DATASET=${DATA_NAME},CLUSTERMETHOD=seurat ../../rnaSeq/code/cluster/cluster_job.sh) 

        echo "seurat job id $cluster2_id"

        start_time=$(date +%M.%S)
        break="false"
        while [ "$break" == "false" ]
        do 
            if [ $(sacct -j ${cluster1_id} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                then
                    if [ $(sacct -j ${cluster2_id} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                        then
                            break="true"
                    fi
            fi
            sleep 30s
            end_time=$(date +%M.%S)
            elapsed=$(echo "scale=0; $end_time - $start_time" | bc)

            echo "Cluster Job Pending time $elapsed "
        done 

        echo "====================== SC3 CONSOLE LOG ======================  " 

        cat "../../rnaSeq/output/cluster/sc3/${DATA_NAME}/sc3-${DATA_NAME}.out"

        echo "====================== Seurat CONSOLE LOG ======================  " 

        cat "../../rnaSeq/output/cluster/seurat/${DATA_NAME}/seurat-${DATA_NAME}.out"
    fi


    if [ "$VISUALIZATION" == "true" ]
    then 
        ###################### VIZUALIZATION JOBS######################

        echo "====================== VIZUALIZATION ANALYSIS ======================  "

        viz1_id=$(sbatch --parsable \
        --job-name=visualize_sc3_tSNE_${DATA_NAME} \
        --output=../../rnaSeq/output/visualize/tSNE/sc3/${DATA_NAME}/tSNE-${DATA_NAME}.out \
        --export=DATASET=${DATA_NAME},CLUSTERMETHOD=sc3,VISUALIZEMETHOD=tSNE ../../rnaSeq/code/visualize/visualize_job.sh)

        echo "sc3 & tSNE job id $viz1_id"

        viz2_id=$(sbatch --parsable \
        --job-name=visualize_seurat_tSNE_${DATA_NAME} \
        --output=../../rnaSeq/output/visualize/tSNE/seurat/${DATA_NAME}/tSNE-${DATA_NAME}.out \
        --export=DATASET=${DATA_NAME},CLUSTERMETHOD=seurat,VISUALIZEMETHOD=tSNE ../../rnaSeq/code/visualize/visualize_job.sh)

        echo "seurat & tSNE job id $viz2_id"

        viz3_id=$(sbatch --parsable \
        --job-name=visualize_sc3_tSNE+PCA_${DATA_NAME} \
        --output=../../rnaSeq/output/visualize/tSNE+PCA/sc3/${DATA_NAME}/tSNE+PCA-${DATA_NAME}.out \
        --export=DATASET=${DATA_NAME},CLUSTERMETHOD=sc3,VISUALIZEMETHOD=tSNE+PCA ../../rnaSeq/code/visualize/visualize_job.sh)

        echo "sc3 & tSNE+PCA job id $viz3_id"

        viz4_id=$(sbatch --parsable \
        --job-name=visualize_seurat_tSNE+PCA_${DATA_NAME} \
        --output=../../rnaSeq/output/visualize/tSNE+PCA/seurat/${DATA_NAME}/tSNE+PCA-${DATA_NAME}.out \
        --export=DATASET=${DATA_NAME},CLUSTERMETHOD=seurat,VISUALIZEMETHOD=tSNE+PCA ../../rnaSeq/code/visualize/visualize_job.sh)

        echo "seurat & tSNE+PCA job id $viz4_id"

        break_viz="False"
        viz_pass="False"

    else 
        break_viz="true"
        viz_pass="true"

    fi

    ###################### COMPARISON JOBS######################

    if [ "$COMPARISON" == "true" ]
    then 
        comp_id1=$(sbatch --parsable \
                    --job-name=AIR+VM+Purity_${DATA_NAME} \
                    --output=../../rnaSeq/output/comparison/AIR+VM+Purity/${DATA_NAME}/AIR+VM+Purity-${DATA_NAME}.out \
                    --export=DATASET=${DATA_NAME},CLUSTERMETHOD=AIR+VM+Purity ../../rnaSeq/code/comparison/comparison_job.sh)
        echo "Job id is $comp_id1 for method: AIR+VM+Purity"


        comp_id2=$(sbatch --parsable \
                    --job-name=chIndex_${DATA_NAME} \
                    --output=../../rnaSeq/output/comparison/chIndex/${DATA_NAME}/chIndex-${DATA_NAME}.out \
                    --export=DATASET=${DATA_NAME},CLUSTERMETHOD=chIndex ../../rnaSeq/code/comparison/comparison_job.sh)
        echo "Job id is $comp_id2 for method: chIndex"

        comp_id3=$(sbatch --parsable \
                    --job-name=cMatrix_${DATA_NAME} \
                    --output=../../rnaSeq/output/comparison/cMatrix/${DATA_NAME}/cMatrix-${DATA_NAME}.out \
                    --export=DATASET=${DATA_NAME},CLUSTERMETHOD=cMatrix ../../rnaSeq/code/comparison/comparison_job.sh)
        echo "Job id is $comp_id3 for method: cMatrix"

        comp_id4=$(sbatch --parsable \
                    --job-name=heatmap_${DATA_NAME} \
                    --output=../../rnaSeq/output/comparison/heatmap/${DATA_NAME}/heatmap-${DATA_NAME}.out \
                    --export=DATASET=${DATA_NAME},CLUSTERMETHOD=heatmap ../../rnaSeq/code/comparison/comparison_job.sh)
        echo "Job id is $comp_id4 for method: heatmap"

        comp_id5=$(sbatch --parsable \
                    --job-name=regGenes_${DATA_NAME} \
                    --output=../../rnaSeq/output/comparison/regGenes/${DATA_NAME}/regGenes-${DATA_NAME}.out \
                    --export=DATASET=${DATA_NAME},CLUSTERMETHOD=regGenes ../../rnaSeq/code/comparison/comparison_job.sh)
        echo "Job id is $comp_id5 for method: regGenes"
        comp_pass="False"
        break_comp="False"

    else 

        comp_pass="true"
        break_comp="true"

    fi


    start_time=$(date +%M.%S)

    while [[ ! ( "$break_viz" == "true" &&"$break_comp" == "true" )  ]]

    do 
    #+++++++++++++ Vizualization Job finish check +++++++++++++#
        if [ "${viz_pass}" == "False" ]
            then ``
                if [ $(sacct -j ${viz1_id} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                    then
                        if [ $(sacct -j ${viz2_id} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                            then
                                if [ $(sacct -j ${viz3_id} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                                    then
                                        if [ $(sacct -j ${viz4_id} --format=State| tail -n1 | xargs) == "COMPLETED" ]
                                            then 
                                                break_viz="true"
                                                viz_pass="True"
                                                echo "All Visualization jobs are complete"
                                        fi
                                fi
                        fi
                fi
        fi

    #+++++++++++++ Comparison Job finish check +++++++++++++#
        if [ "${comp_pass}" == "False" ]
            then 
                if [ $(sacct -j ${comp_id1} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                    then
                        if [ $(sacct -j ${comp_id2} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                            then
                                if [ $(sacct -j ${comp_id3} --format=State| tail -n1 | xargs) == "COMPLETED" ] 
                                    then
                                        if [ $(sacct -j ${comp_id4} --format=State| tail -n1 | xargs) == "COMPLETED" ]
                                            then 
                                                if [ $(sacct -j ${comp_id5} --format=State| tail -n1 | xargs) == "COMPLETED" ]
                                                    then 
                                                        break_comp="true"
                                                        comp_pass="True"
                                                        echo "All Comparison jobs are complete"
                                                fi
                                        fi
                                fi
                        fi
                fi
        fi
        sleep 30s
        end_time=$(date +%M.%S)
        elapsed=$(echo "scale=0; $end_time - $start_time" | bc)

        echo "Job Pending time $elapsed "
    done 
    if [ "$COMPARISON" == "true" ]
    then 
        echo "====================== Comparison CONSOLE LOG ======================  " 
        cat ../../rnaSeq/output/comparison/AIR+VM+Purity/${DATA_NAME}/AIR+VM+Purity-${DATA_NAME}.out
        cat ../../rnaSeq/output/comparison/chIndex/${DATA_NAME}/chIndex-${DATA_NAME}.out
        cat ../../rnaSeq/output/comparison/cMatrix/${DATA_NAME}/cMatrix-${DATA_NAME}.out
        cat ../../rnaSeq/output/comparison/heatmap/${DATA_NAME}/heatmap-${DATA_NAME}.out
        cat ../../rnaSeq/output/comparison/regGenes/${DATA_NAME}/regGenes-${DATA_NAME}.out
    fi
    if [ "$VISUALIZATION" == "true" ]
    then    
        echo "====================== Visualization CONSOLE LOG ======================  " 
        cat ../../rnaSeq/output/visualize/tSNE+PCA/seurat/${DATA_NAME}/tSNE+PCA-${DATA_NAME}.out
        cat ../../rnaSeq/output/visualize/tSNE+PCA/sc3/${DATA_NAME}/tSNE+PCA-${DATA_NAME}.out
        cat ../../rnaSeq/output/visualize/tSNE/seurat/${DATA_NAME}/tSNE-${DATA_NAME}.out
        cat ../../rnaSeq/output/visualize/tSNE/seurat/${DATA_NAME}/tSNE-${DATA_NAME}.out
    fi


    git add --all
    git commit -m "Ran PREPROCESS : $PREPROCESS 
    > Ran CLUSTER : $CLUSTER 
    > Ran VISUALIZATION : $VISUALIZATION 
    > Ran COMPARISON : $COMPARISON"
    git push 


