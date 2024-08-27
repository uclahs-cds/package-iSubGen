# test.R
# this is a test script updated based on the README.md to test all the functions in the package

library(iSubGen);
library(tensorflow);
library(reticulate);
library(keras3);

## testing load.molecular.aberration.data() ########################################
molecular.data <- list();
for(i in c('cna','methy','snv')) {
    molecular.data[[i]] <- load.molecular.aberration.data(
        system.file('exdata',paste0(i,'_profiles.txt'), package='iSubGen'),
        patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
        );
}

## testing calculate.scaling() #####################################################
# scale the mrna and mirna data
scaling.factors <- calculate.scaling(
	list(
		mirna = molecular.data$mirna,
		mrna_tc_together = molecular.data$mrna_tc_together,
		mrna_tac_together = molecular.data$mrna_tac_together
		)
	);

## testing write.scaling.factors() #################################################
# write the scaling factors to a file
write.scaling.factors(
	scaling.factors,
	paste0(pipe.vars$pipeline.dir,'input_data/scaling_factors/')
	);

## testing apply.scaling() #########################################################
scaled.matrices <- apply.scaling(
	list(
		mirna = molecular.data$mirna,
		mrna_tc_together = molecular.data$mrna_tc_together,
		mrna_tac_together = molecular.data$mrna_tac_together
		),
	scaling.factors
	);
molecular.data$mrna_tc_together <- scaled.matrices$mrna_tc_together;
molecular.data$mrna_tac_together <- scaled.matrices$mrna_tac_together;
molecular.data$mirna <- scaled.matrices$mirna;

## testing create.autoencoder() ####################################################
# Create a list to store the autoencoders
autoencoders <- list();

# Create and train an autoencoder using CNA data
autoencoders[['cna']] <- create.autoencoder(
    data.type = 'cna',
    data.matrix = molecular.data$cna,
    encoder.layers.node.nums = c(20,2)
    )$autoencoder;
# take a look at the layers/number of nodes in the autoencoder
str(autoencoders$cna);

# Create and train an autoencoder using methylation data
autoencoders[['methy']] <- create.autoencoder(
    data.type = 'methy',
    data.matrix = molecular.data$methy,
    encoder.layers.node.nums = c(15,1)
    )$autoencoder;
# take a look at the layers/number of nodes in the autoencoder
str(autoencoders$methy);

# Create and train an autoencoder using coding SNV data
autoencoders[['snv']] <- create.autoencoder(
    data.type = 'snv',
    data.matrix = molecular.data$snv,
    encoder.layers.node.nums = c(15,1)
    )$autoencoder;

## test create.autoencoder.irf.matrix() ############################################
# Get the independent reduced features from the autoencoders
irf.matrix <- create.autoencoder.irf.matrix(
    data.types = names(molecular.data),
    data.matrices = molecular.data,
    autoencoders = autoencoders
    );

## test calculate.integrative.similarity.matrix() ##################################
# Calculate a similarity matrix using correlations
similarity.matrix <- calculate.integrative.similarity.matrix(
    data.types = names(molecular.data),
    data.matrices = molecular.data,
    dist.metrics = list(
        cna = 'euclidean',
        snv = 'euclidean',
        methy = 'euclidean'
        )
    );

## test calculate.cis.matrix() #####################################################
cis.matrix <- calculate.cis.matrix(
    data.types = names(molecular.data),
    data.matrices = molecular.data,
    dist.metrics = list(
        cna = 'euclidean',
        snv = 'euclidean',
        methy = 'euclidean'
        )
    );

## test combine.integrative.features() ##############################################
# Combine IRFs and CISs into one matrix
integrative.features.matrix <- combine.integrative.features(
    irf.matrix,
    cis.matrix
    )$integrative.feature.matrix;

## test cluster.patients() #########################################################
# Perform consensus clustering to get integrative subtypes
subtyping.results <- cluster.patients(
    data.matrix = integrative.features.matrix,
    distance.metric = 'euclidean',
    parent.output.dir = './',
    new.result.dir = 'vignette_subtypes',
    max.num.subtypes = 5,
    clustering.reps = 50,
    consensus.cluster.write.table = FALSE
    );
