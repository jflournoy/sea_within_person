#!/bin/bash
3dUndump -ijk -datum 'float' -prefix TRR_meta_cont_est.nii.gz	-master Sch400_HO_combined.nii.gz TRR_meta_cont_est.txt
3dUndump -ijk -datum 'float' -prefix TRR_meta_cont_ll.nii.gz	-master Sch400_HO_combined.nii.gz TRR_meta_cont_ll.txt
3dUndump -ijk -datum 'float' -prefix TRR_meta_cont_uu.nii.gz	-master Sch400_HO_combined.nii.gz TRR_meta_cont_uu.txt
3dUndump -ijk -datum 'float' -prefix TRR_meta_avg_est.nii.gz 	-master Sch400_HO_combined.nii.gz TRR_meta_avg_est.txt
3dUndump -ijk -datum 'float' -prefix TRR_meta_avg_ll.nii.gz 	-master Sch400_HO_combined.nii.gz TRR_meta_avg_ll.txt
3dUndump -ijk -datum 'float' -prefix TRR_meta_avg_uu.nii.gz 	-master Sch400_HO_combined.nii.gz TRR_meta_avg_uu.txt

3dUndump -ijk -datum 'float' -prefix oldICC_meta_est.nii.gz	-master Sch400_HO_combined.nii.gz oldICC_meta_est.txt 
3dUndump -ijk -datum 'float' -prefix oldICC_meta_ll.nii.gz	-master Sch400_HO_combined.nii.gz oldICC_meta_ll.txt
3dUndump -ijk -datum 'float' -prefix oldICC_meta_uu.nii.gz	-master Sch400_HO_combined.nii.gz oldICC_meta_uu.txt
3dUndump -ijk -datum 'float' -prefix oldICC_all_est.nii.gz 	-master Sch400_HO_combined.nii.gz oldICC_all_est.txt
3dUndump -ijk -datum 'float' -prefix oldICC_all_ll.nii.gz 	-master Sch400_HO_combined.nii.gz oldICC_all_ll.txt
3dUndump -ijk -datum 'float' -prefix oldICC_all_uu.nii.gz 	-master Sch400_HO_combined.nii.gz oldICC_all_uu.txt
