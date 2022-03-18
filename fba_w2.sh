Cross-sectional FBA on only Wave 2 people for the follow-up vs. clinically longitudinal Wave1-to-Wave2 models (N=115)
#https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/st_fibre_density_cross-section.html

# set up environment
module load mrtrix
module load ANTs
module load anaconda3/5.0.0
source activate trix # see: dwifslpreproc_conda_HPC.sh

# define user vars
W2_DIR=/project/3022043.01/NeuroIMAGE2/niftibox
FBA_DIR=/project/3022028.01/fixel_based_analysis/fba_w2
LONG_DIR=/project/3022028.01/fixel_based_analysis/fba_long_w1_w2
SUBJECT_LIST=`cat $FBA_DIR/subject_list_w2.txt`
TEMPLATE_LIST=`cat $FBA_DIR/template_list_w2.txt`

cd $FBA_DIR

for IID in $SUBJECT_LIST;do

# create output directories for wave 2
mkdir -p $FBA_DIR/${IID}

# NeuroIMAGE 2 data have already been preprocessed a while ago, so the data need only to be format converted from nifti to mif
## Convert pre-existing preprocessed diffusion data into .mif
mrconvert $W2_DIR/${IID}/diffusion/FDT_Data/data.nii.gz $FBA_DIR/${IID}/dwi.mif -fslgrad $W2_DIR/${IID}/diffusion/FDT_Data/bvecs $W2_DIR/${IID}/diffusion/FDT_Data/bvals -force
done

# 3. Estimate a temporary brain mask
for IID in $SUBJECT_LIST;do
echo "dwi2mask $FBA_DIR/${IID}/dwi.mif $FBA_DIR/${IID}/dwi_temp_mask.mif -force" |qsub -N dwi2mask -V -l walltime=00:05:00,mem=10gb
done

# 4. Bias field correction
for IID in $SUBJECT_LIST;do
echo "dwibiascorrect ants $FBA_DIR/${IID}/dwi.mif $FBA_DIR/${IID}/dwi_unbiased.mif -force" | qsub -N dwibiascorrect -V -l walltime=00:10:00,mem=10gb
done

# 5. Global intensity normalisation across subjects
##  create directories to store all the input and output images
mkdir -p $FBA_DIR/dwinormalise/dwi_input
mkdir $FBA_DIR/dwinormalise/mask_input

for IID in $SUBJECT_LIST;do
ln -sr $FBA_DIR/${IID}/dwi_unbiased.mif $FBA_DIR/dwinormalise/dwi_input/${IID}.mif
ln -sr $FBA_DIR/${IID}/dwi_temp_mask.mif $FBA_DIR/dwinormalise/mask_input/${IID}.mif
done
## perform group DWI intensity normalisation in cluster
echo "dwinormalise group $FBA_DIR/dwinormalise/dwi_input/ $FBA_DIR/dwinormalise/mask_input/ $FBA_DIR/dwinormalise/dwi_output/ $FBA_DIR/dwinormalise/fa_template.mif $FBA_DIR/dwinormalise/fa_template_wm_mask.mif -force" | qsub -N dwinormalise -V -l walltime=06:00:00,mem=20gb

#### edit non-brain out of fa_template_wm_mask 
## link the output files back to the subject directories
for IID in $SUBJECT_LIST;do
ln -sr -f $FBA_DIR/dwinormalise/dwi_output/${IID}.mif $FBA_DIR/${IID}/dwi_normalised.mif
done

#6. Computing (average) tissue response functions
for IID in $SUBJECT_LIST;do
echo "dwi2response tournier $FBA_DIR/${IID}/dwi_normalised.mif $FBA_DIR/${IID}/response.txt -force" | qsub -N dwi2response -V -l walltime=00:10:00,mem=5gb
done
### OR multi-tissue response
for IID in $SUBJECT_LIST;do
echo "dwi2response dhollander $FBA_DIR/${IID}/dwi.mif $FBA_DIR/${IID}/response_wm.txt $FBA_DIR/${IID}/response_gm.txt $FBA_DIR/${IID}/response_csf.txt -force" | qsub -N dwi2response -V -l walltime=00:10:00,mem=5gb
done

## average the response functions obtained from all subjects
responsemean $FBA_DIR/*/response.txt $FBA_DIR/group_average_response.txt -force
### OR 2-tissue response
responsemean */response_wm.txt group_average_response_wm.txt -force
responsemean */response_csf.txt group_average_response_csf.txt -force

# 7. Upsampling DW images 
for IID in $SUBJECT_LIST;do
echo "mrgrid $FBA_DIR/${IID}/dwi_normalised.mif regrid -vox 1.25 $FBA_DIR/${IID}/dwi_normalised_upsampled.mif -force" | qsub -N mrgrid -V -l walltime=00:05:00,mem=10gb
done

### OR 2-tissue response
for IID in $SUBJECT_LIST;do
echo "mrgrid $FBA_DIR/${IID}/dwi_unbiased.mif regrid -vox 1.25 $FBA_DIR/${IID}/dwi_unbiased_upsampled.mif -force" | qsub -N mrgrid -V -l walltime=00:05:00,mem=10gb
done

# 8. Compute upsampled whole brain mask images
for IID in $SUBJECT_LIST;do
echo "dwi2mask $FBA_DIR/${IID}/dwi_normalised_upsampled.mif $FBA_DIR/${IID}/dwi_mask_upsampled.mif -force" | qsub -N dwi2mask -V -l walltime=00:05:00,mem=5gb
done

# 9. Fibre Orientation Distribution estimation (spherical deconvolution)
for IID in $SUBJECT_LIST;do
echo "dwiextract $FBA_DIR/${IID}/dwi_normalised_upsampled.mif $FBA_DIR/${IID}/dwi_normalised_upsampled_dwiextract.mif -force
dwi2fod -mask $FBA_DIR/${IID}/dwi_mask_upsampled.mif msmt_csd $FBA_DIR/${IID}/dwi_normalised_upsampled_dwiextract.mif $FBA_DIR/group_average_response.txt $FBA_DIR/${IID}/wmfod.mif -force" | qsub -N dwiextract2fod -V -l walltime=02:00:00,mem=5gb
done
## OR 2-tissue multi-tissue CSD
# https://community.mrtrix.org/t/fba-single-shell-high-gm-fod-values/1746/2
for IID in $SUBJECT_LIST;do
echo "dwi2fod msmt_csd $FBA_DIR/${IID}/dwi_unbiased_upsampled.mif $FBA_DIR/group_average_response_wm.txt $FBA_DIR/${IID}/wmfod.mif $FBA_DIR/group_average_response_csf.txt $FBA_DIR/${IID}/csf.mif -mask $FBA_DIR/${IID}/dwi_mask_upsampled.mif -force" | qsub -N dwi2fod -V -l walltime=02:00:00,mem=5gb
done

## OR 2-tissue Joint bias field correction and intensity normalisation
for IID in $SUBJECT_LIST;do
echo "mtnormalise $FBA_DIR/${IID}/wmfod.mif $FBA_DIR/${IID}/wmfod_norm.mif $FBA_DIR/${IID}/csf.mif $FBA_DIR/${IID}/csf_norm.mif -mask $FBA_DIR/${IID}/dwi_mask_upsampled.mif -force" | qsub -N mtnormalise -V -l walltime=02:00:00,mem=5gb
done

# 10. Generate a study-specific unbiased FOD template
mkdir -p $FBA_DIR/template/fod_input
mkdir $FBA_DIR/template/mask_input
## Symbolic link FOD images (and masks) into a single input folder
for IID in $SUBJECT_LIST;do
ln -sr $FBA_DIR/${IID}/wmfod_norm.mif $FBA_DIR/template/fod_input/${IID}.mif
ln -sr $FBA_DIR/${IID}/dwi_mask_upsampled.mif $FBA_DIR/template/mask_input/${IID}.mif
done
## Run the template building script 
echo "population_template $FBA_DIR/template/fod_input -mask_dir $FBA_DIR/template/mask_input $FBA_DIR/template/wmfod_template.mif -voxel_size 1.25 -scratch $FBA_DIR/template/ -continue $FBA_DIR/template/population_template-tmp-Y541N3/ warps_05/518-1377.mif -nocleanup -force > $FBA_DIR/template/population_template.log 2>&1" | qsub -N wmfod_template -V -l nodes=1:ppn=9,walltime=48:00:00,mem=256gb
### process will time out of cluster, so need to use -continue option the next time around, using the last successfully created file (last shared file between two warp dirs)
-continue $FBA_DIR/template/population_template-tmp-Y541N3/ warps_05/518-1377.mif

# 11. Register all subject FOD images to the FOD template
for IID in $SUBJECT_LIST;do 
echo "mrregister $FBA_DIR/${IID}/wmfod.mif -mask1 $FBA_DIR/${IID}/dwi_mask_upsampled.mif $FBA_DIR/template/wmfod_template.mif -nl_warp $FBA_DIR/${IID}/subject2template_warp.mif $FBA_DIR/${IID}/template2subject_warp.mif -force" | qsub -N mrregister -V -l walltime=02:00:00,mem=12gb
done

# 12. Compute the template mask (intersection of all subject masks in template space)
for IID in $SUBJECT_LIST;do 
echo "mrtransform $FBA_DIR/${IID}/dwi_mask_upsampled.mif -warp $FBA_DIR/${IID}/subject2template_warp.mif -interp nearest -datatype bit $FBA_DIR/${IID}/dwi_mask_in_template_space.mif -force" | qsub -N mrtransform -V -l walltime=00:05:00,mem=2gb
done
#Compute the template mask as the intersection of all warped masks
mrmath $FBA_DIR/*/dwi_mask_in_template_space.mif min $FBA_DIR/template/template_mask.mif -datatype bit -force

# 13. Compute a white matter template analysis fixel mask
fod2fixel -mask $FBA_DIR/template/template_mask.mif -fmls_peak_value 0.10 $FBA_DIR/template/wmfod_template.mif $FBA_DIR/template/fixel_mask

# 14. Warp FOD images to template space
for IID in $SUBJECT_LIST;do 
echo "mrtransform $FBA_DIR/${IID}/wmfod.mif -warp $FBA_DIR/${IID}/subject2template_warp.mif -reorient_fod no $FBA_DIR/${IID}/fod_in_template_space_NOT_REORIENTED.mif -force" | qsub -N mrtransform -V -l walltime=00:05:00,mem=2gb
done

# 15. Segment FOD images to estimate fixels and their apparent fibre density (FD)
for IID in $SUBJECT_LIST;do 
echo "fod2fixel -mask $FBA_DIR/template/template_mask.mif $FBA_DIR/${IID}/fod_in_template_space_NOT_REORIENTED.mif $FBA_DIR/${IID}/fixel_in_template_space_NOT_REORIENTED -afd fd.mif -force" | qsub -N fod2fixel -V -l walltime=01:00:00,mem=15gb
done

# 16. Reorient fixels
for IID in $SUBJECT_LIST;do 
echo "fixelreorient $FBA_DIR/${IID}/fixel_in_template_space_NOT_REORIENTED $FBA_DIR/${IID}/subject2template_warp.mif $FBA_DIR/${IID}/fixel_in_template_space -force" | qsub -N fixelreorient -V -l walltime=00:05:00,mem=15gb
done
rm -r */fixel_in_template_space_NOT_REORIENTED

# 17. Assign subject fixels to template fixels
for IID in $SUBJECT_LIST;do 
fixelcorrespondence $FBA_DIR/${IID}/fixel_in_template_space/fd.mif $FBA_DIR/template/fixel_mask $FBA_DIR/template/fd ${IID}.mif -force 
done

# 18. Compute the fibre cross-section (FC) metric
for IID in $SUBJECT_LIST;do 
warp2metric $FBA_DIR/${IID}/subject2template_warp.mif -fc $FBA_DIR/template/fixel_mask $FBA_DIR/template/fc ${IID}.mif -force
done

mkdir $FBA_DIR/template/log_fc
cp $FBA_DIR/template/fc/index.mif $FBA_DIR/template/fc/directions.mif $FBA_DIR/template/log_fc
for IID in $SUBJECT_LIST;do 
mrcalc $FBA_DIR/template/fc/${IID}.mif -log $FBA_DIR/template/log_fc/${IID}.mif -force
done

# 19. Compute a combined measure of fibre density and cross-section (FDC)
mkdir $FBA_DIR/template/fdc
cp $FBA_DIR/template/fc/index.mif $FBA_DIR/template/fdc
cp $FBA_DIR/template/fc/directions.mif $FBA_DIR/template/fdc
for IID in $SUBJECT_LIST;do 
mrcalc $FBA_DIR/template/fd/${IID}.mif $FBA_DIR/template/fc/${IID}.mif -mult $FBA_DIR/template/fdc/${IID}.mif -force
done
### mrcalc: [WARNING] header transformations of input images do not match

# 20. Perform whole-brain fibre tractography on the FOD template
cd $FBA_DIR/template
tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_image template_mask.mif -mask template_mask.mif -select 20000000 -cutoff 0.15 tracks_20_million.tck -force

## The appropriate FOD amplitude cutoff for FOD template tractography can vary considerably between different datasets, as well as different versions of MRtrix3 due to historical software bugs. While the value of 0.10 is suggested as a reasonable value for single-tissue data, it may be beneficial to first generate a smaller number of streamlines (e.g. 100,000) using this value, and visually confirm that the generated streamlines exhibit an appropriate extent of propagation at the ends of white matter pathways, before committing to generation of the dense tractogram.

# 21. Reduce biases in tractogram densities
echo "tcksift $FBA_DIR/template/tracks_20_million.tck $FBA_DIR/template/wmfod_template.mif $FBA_DIR/template/tracks_2_million_sift.tck -term_number 2000000 -force" | qsub -N tcksift -V -l walltime=05:00:00,mem=40gb

# 22. Generate fixel-fixel connectivity matrix
echo "fixelconnectivity $FBA_DIR/template/fixel_mask/ $FBA_DIR/template/tracks_2_million_sift.tck $FBA_DIR/template/matrix/ -force" | qsub -N fixelconnect -V -l walltime=00:15:00,mem=20gb

# 23. Smooth fixel data using fixel-fixel connectivity
for F in fd log_fc fdc;do
echo "fixelfilter $FBA_DIR/template/${F} smooth $FBA_DIR/template/${F}_smooth -matrix $FBA_DIR/template/matrix/ -force" | qsub -N fixelfilter_${F} -V -l walltime=00:45:00,mem=8gb
done

# Run tractseg.sh to generate ROI masks for lCST and lSLF

# rm excluded subjects (bc no CPRS w2 data) from fixel dirs, design matrix, and files.txt
#cd /project/3022028.01/fixel_based_analysis/fba_w2/template
#for dir in fd_smooth fdc_smooth log_fc_smooth;do
#for IID in 98-114-10059 98-128-10117 98-131-10131 98-139-10166 98-171-10307 98-174-10320 98-181-10452 98-228-20134 98-308-30034 98-308-30035 98-333-30142 98-341-30175 98-350-30214 98-350-30215 98-527-1414 98-539-1467;do
#rm -r ${IID}
#done
#done

# 24. Perform statistical analysis of FD, FC, and FDC
# write scripts for fixelcfestats
TEMPLATE_DIR=/project/3022028.01/fixel_based_analysis/fba_w2/template
OUT_DIR=$LONG_DIR/fixelcfestats
cd $OUT_DIR

# create files.txt

# total symptom score change
for F in fdc_smooth fd_smooth log_fc_smooth;do
for D in design_total_demean;do
for T in CST_left SLF_I_left SLF_II_left SLF_III_left;do
mkdir -p $OUT_DIR/fixelcfestats_out/21jul2021
command_CFE="fixelcfestats -mask /project/3022028.01/fixel_based_analysis/tractseg_w2/tck2fixel_out/${T}.mif -permutations $OUT_DIR/palm/permutations_set_5000.csv $TEMPLATE_DIR/${F} $TEMPLATE_DIR/files.txt $OUT_DIR/${D}.csv $OUT_DIR/contrast.txt $TEMPLATE_DIR/matrix $OUT_DIR/fixelcfestats_out/21jul2021/${F}_${D}_${T}/ -force" 
SCRIPT="$OUT_DIR/fixelcfestats_out/21jul2021/script_${F}_${D}_${T}"
echo "$command_CFE" > $SCRIPT
chmod a+rwx $SCRIPT
done
done
done

# hyp-imp score change 
for F in fdc_smooth fd_smooth log_fc_smooth;do
for D in design_HI_demean;do
for T in CST_left;do
mkdir -p $OUT_DIR/fixelcfestats_out/21jul2021
command_CFE="fixelcfestats -mask /project/3022028.01/fixel_based_analysis/tractseg_w2/tck2fixel_out/${T}.mif -permutations $OUT_DIR/palm/permutations_set_5000.csv $TEMPLATE_DIR/${F} $TEMPLATE_DIR/files.txt $OUT_DIR/${D}.csv $OUT_DIR/contrast.txt $TEMPLATE_DIR/matrix $OUT_DIR/fixelcfestats_out/21jul2021/${F}_${D}_${T}/ -force" 
SCRIPT="$OUT_DIR/fixelcfestats_out/21jul2021/script_${F}_${D}_${T}"
echo "$command_CFE" > $SCRIPT
chmod a+rwx $SCRIPT
done
done
done

# submit each script to cluster
cd $OUT_DIR/fixelcfestats_out/21jul2021
for script in script_fdc_smooth_design_HI_demean_CST_left script_fd_smooth_design_total_demean_SLF_II_left script_fdc_smooth_design_total_demean_CST_left      script_fd_smooth_design_total_demean_SLF_I_left script_fdc_smooth_design_total_demean_SLF_III_left  script_log_fc_smooth_design_HI_demean_CST_left script_fdc_smooth_design_total_demean_SLF_II_left   script_log_fc_smooth_design_total_demean_CST_left script_fdc_smooth_design_total_demean_SLF_I_left    script_log_fc_smooth_design_total_demean_SLF_III_left script_fd_smooth_design_HI_demean_CST_left          script_log_fc_smooth_design_total_demean_SLF_II_left script_fd_smooth_design_total_demean_CST_left       script_log_fc_smooth_design_total_demean_SLF_I_left script_fd_smooth_design_total_demean_SLF_III_left;do
qsub -N fixelcfestats -l walltime=24:00:00,mem=64gb ${script}
done

# view results as streamlines
fixel2tsf fixelcfestats_out/21jul2021/fd_smooth_design_HI_demean_CST_left/fwe_1mpvalue_t2.mif $TEMPLATE_DIR/tracks_200k_sift.tck fd_fwe_pvalue.tsf

fixel2tsf fixelcfestats_out/21jul2021/fd_smooth_design_HI_demean_CST_left/std_effect_t2.mif $TEMPLATE_DIR/tracks_200k_sift.tck fd_std_effect.tsf




