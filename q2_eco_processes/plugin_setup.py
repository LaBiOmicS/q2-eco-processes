import qiime2.plugin
from qiime2.plugin import Plugin, Metadata, Str, Int
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.tree import Phylogeny, Rooted
from q2_types.distance_matrix import DistanceMatrix

import q2_eco_processes
from q2_eco_processes._processes import calculate_processes, calculate_matrices
from q2_eco_processes._visualizer import summarize_processes
from q2_eco_processes._signal import check_phylogenetic_signal

plugin = Plugin(
    name='eco-processes',
    version=q2_eco_processes.__version__,
    website='https://github.com/LaBiOmics/q2-eco-processes',
    package='q2_eco_processes',
    user_support_text='https://github.com/LaBiOmics/q2-eco-processes/issues',
    citation_text=(
        "Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ, et al. (2013) Quantifying community assembly processes and identifying features that impose them. The ISME Journal 7(11):2069-2079."
    ),
    description=(
        'QIIME 2 plugin for quantifying ecological assembly processes (Homogeneous Selection, '
        'Variable Selection, Dispersal Limitation, Homogenizing Dispersal, and Undominated Drift) '
        'using phylogenetic null models across experimental groups.'
    ),
    short_description='Microbiome ecological assembly processes plugin.'
)

# Register Method Action: calculate-processes
plugin.methods.register_function(
    function=calculate_processes,
    inputs={
        'table': FeatureTable[Frequency],
        'tree': Phylogeny[Rooted]
    },
    parameters={
        'metadata': Metadata,
        'column': Str,
        'permutations': Int,
        'min_frequency': Int,
        'n_jobs': Int
    },
    outputs=[
        ('process_table', FeatureTable[Frequency])
    ],
    input_descriptions={
        'table': 'ASV feature table (FeatureTable[Frequency]).',
        'tree': 'Rooted phylogenetic tree (Phylogeny[Rooted]).'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata for experimental group comparisons (Control vs Treated).',
        'column': 'Metadata column name defining experimental groups.',
        'permutations': 'Number of null model randomization permutations (default: 999).',
        'min_frequency': 'Filter out low-frequency ASVs with total count below threshold (default: 0).',
        'n_jobs': 'Number of CPU threads to use for parallel permutations (default: 1).'
    },
    output_descriptions={
        'process_table': 'Ecological assembly process frequency table per sample or group (FeatureTable[Frequency]).'
    },
    name='Ecological Assembly Processes Quantification',
    description='Quantifies the relative influence of selection, dispersal, and drift using Stegen null models.'
)

# Register Method Action: calculate-matrices
plugin.methods.register_function(
    function=calculate_matrices,
    inputs={
        'table': FeatureTable[Frequency],
        'tree': Phylogeny[Rooted]
    },
    parameters={
        'permutations': Int,
        'min_frequency': Int,
        'n_jobs': Int
    },
    outputs=[
        ('beta_nti', DistanceMatrix),
        ('rc_bray', DistanceMatrix)
    ],
    input_descriptions={
        'table': 'ASV feature table (FeatureTable[Frequency]).',
        'tree': 'Rooted phylogenetic tree (Phylogeny[Rooted]).'
    },
    parameter_descriptions={
        'permutations': 'Number of null model randomization permutations (default: 999).',
        'min_frequency': 'Filter out low-frequency ASVs with total count below threshold (default: 0).',
        'n_jobs': 'Number of CPU threads to use for parallel permutations (default: 1).'
    },
    output_descriptions={
        'beta_nti': 'Pairwise BetaNTI score matrix (DistanceMatrix).',
        'rc_bray': 'Pairwise Raup-Crick (RCbray) score matrix (DistanceMatrix).'
    },
    name='Ecological Process Pairwise Distance Matrices',
    description='Calculates pairwise BetaNTI and Raup-Crick (RCbray) matrices as DistanceMatrix artifacts.'
)

# Register Visualizer Action: summarize-processes
plugin.visualizers.register_function(
    function=summarize_processes,
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'metadata': Metadata
    },
    input_descriptions={
        'table': 'Process table produced by calculate-processes.'
    },
    parameter_descriptions={
        'metadata': 'Optional sample metadata for sample grouping and metadata coloring.'
    },
    name='Ecological Assembly Process Summary Visualization',
    description='Generates an interactive HTML summary dashboard (.qzv) summarizing ecological assembly processes.'
)

# Register Visualizer Action: check-phylogenetic-signal
plugin.visualizers.register_function(
    function=check_phylogenetic_signal,
    inputs={
        'table': FeatureTable[Frequency],
        'tree': Phylogeny[Rooted]
    },
    parameters={
        'metadata': Metadata,
        'column': Str,
        'permutations': Int
    },
    input_descriptions={
        'table': 'ASV feature table (FeatureTable[Frequency]).',
        'tree': 'Rooted phylogenetic tree (Phylogeny[Rooted]).'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata containing numeric environmental variables.',
        'column': 'Numeric metadata column name to evaluate niche distance (e.g. pH, Temperature).',
        'permutations': 'Number of Mantel test permutations (default: 999).'
    },
    name='Phylogenetic Signal Diagnostic Visualization',
    description='Evaluates phylogenetic niche conservatism by testing correlation between niche distance and phylogenetic distance.'
)
