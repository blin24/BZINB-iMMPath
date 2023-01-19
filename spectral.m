filename = "cut_weights.txt";
delimiterIn = ' ';
headerlinesIn = 0;
T = importdata(filename,delimiterIn,headerlinesIn);


SPECTRAL_HOME = '/proj/epi/CVDGeneNas/bmlin/microbiome/cut/SpectraLib'
ASYMSPECTRAL_HOME = '/proj/epi/CVDGeneNas/bmlin/microbiome/cut/SpectraLib_A'
scriptname='/proj/epi/CVDGeneNas/bmlin/microbiome/cut/SpectraLib_A/install_SpectraLib_A.m'
run(scriptname)




filename = "affinity_bzinb_microb_ecc0.txt";
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);

filename = "affinity_bzinb_microb_ecc1.txt";
delimiterIn = ' ';
headerlinesIn = 0;
A2 = importdata(filename,delimiterIn,headerlinesIn);


asig=cluster_wcut(A,T,6)
asig2=cluster_wcut(A2,T,6)



T = table(asig', asig2', 'VariableNames', { 'ECC0', 'ECC1'} )

writetable(T, 'out.txt')


