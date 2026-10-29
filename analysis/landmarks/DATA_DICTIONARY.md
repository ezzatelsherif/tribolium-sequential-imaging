# Landmark data columns

Blank numeric fields indicate unavailable measurements or intervals. Spatial positions are fractions of whole embryo length, increasing from anterior to posterior. Temporal positions and differences are ordered-substage coordinates, not elapsed time.

## S1 Data: Landmark inventory

| Column | Meaning |
| --- | --- |
| `comparison_id` | Unique identifier for a gene/event pair, spatial stage and landmark. |
| `contrast_id` | Temporal event-pair/landmark identifier. The same temporal estimate can recur at several spatial stages. |
| `group` | `pair-rule`, `gap-gap`, or `gap-pair-rule`. |
| `landmark` | `peak`, `rise`, or `fall`. |
| `stage` | Spatial stage: `Stage3_3`, `Stage4_1`, `Stage4_3`, or `Stage5_2`. |
| `gene_a`, `gene_b` | The two genes in the canonical A/B order. |
| `event_a`, `event_b` | Assigned expression episodes or pair-rule cycles. `gt_episode` is the shared outer episode used for gt borders. |
| `cohort_a`, `cohort_b` | Temporal cohort identities, preserving target-matched eve channels. |
| `pair_label` | Human-readable event-pair label. |
| `time_a`, `time_b` | Observed temporal landmark estimates. |
| `dt` | Canonical temporal difference, B minus A. Positive means B occurs later. |
| `ci_low`, `ci_high` | Pointwise 95% percentile interval for the direct B-minus-A temporal difference. |
| `time_a_ci_high`, `time_b_ci_high` | Upper marginal 95% temporal limits used by the common availability criterion. |
| `matched_temporal_channels` | True when both temporal measurements belong to the same jointly resampled cohort. |
| `ds_emb1`, `ds_emb2`, `ds_emb3` | B-minus-A spatial differences in each embryo. Positive means B is posterior to A. |
| `position_a_emb1` through `position_a_emb3` | Spatial positions of landmark A in the three embryos. |
| `position_b_emb1` through `position_b_emb3` | Spatial positions of landmark B in the three embryos. |
| `eligible` | True only if the entry passes all temporal, spatial and event-identity criteria. |
| `exclusion_reasons` | Semicolon-separated reasons for excluding an entry. All excluded candidates are retained. |
| `flags` | Additional source feature and event-identity annotations. |
| `recovery` | Fraction of 20,000 bootstrap draws with a valid direct temporal contrast. |
| `outcome` | Original classification before separating consistent spatial ordering with unresolved timing. |
| `span_outcome` | Alternative peak classification using daughter-peak spans for split eve, rather than only their midpoint. |
| `orientation` | Sign of the observed temporal difference, retained from the source table; not used to change classifications between bootstrap draws. |
| `peak_definition` | `single_peak`, `split_eve_midpoint`, or `outer_episode_halfheight`. |
| `temporal_source_a`, `temporal_source_b` | Keys of the corresponding arrays in the saved temporal feature NPZ. |
| `summary_panel` | `all_landmarks` for eligible comparisons in Fig. 7B, or `excluded`. The original `landmark` field retains peak, rise or fall. |
| `summary_category` | One of the five plotted categories for eligible entries, or `excluded`. |
| `summary_category_label` | Human-readable category label, including `Unresolved in time, separated in space`. |

A positive temporal and positive spatial difference preserve order. A negative temporal and negative spatial difference also preserve order. Reversal requires opposite signs, a temporal interval excluding zero, and agreement in all three spatial embryos.

## S1 Data: Summary source

| Column | Meaning |
| --- | --- |
| `panel` | `all_landmarks` for the pooled figure summary, `peak` for the peak stratum, or `half_height` for the combined rising/falling stratum. |
| `group` | Gene-pair group using the same labels as the Landmark inventory sheet. |
| `category` | One of the five plotted categories. Every category is retained, including zero counts. |
| `count` | Number of eligible comparison entries in this category. |
| `total` | Total eligible entries for that panel and gene-pair group. |
| `percentage` | `100 * count / total`. Values are percentages on a 0-100 scale. |
| `summary_scope` | `pooled_landmarks`, `peaks`, or `half_heights`. |
| `figure_panel` | `B` for the 15 plotted pooled rows; blank for the 30 supplementary stratified rows. |
| `category_label` | Human-readable outcome label. |

The table contains 45 rows: 15 pooled rows for Figure 7B and 30 retained stratified rows. Each summary has three gene-pair groups and five outcome categories, including zeros. The pooled rows combine eligible peak, rising half-height and falling half-height comparisons. The stratified rows show the same observations as peaks versus combined rising/falling half-heights for supplementary reporting. Do not sum pooled and stratified rows together. Percentages are calculated from the corresponding counts and denominators, not by averaging separate percentages.

For `spatial_consistent_temporal_unresolved`, the requested label is **Unresolved in time, separated in space**. This means the 95% temporal difference interval includes zero and all three spatial differences have the same nonzero sign. It does not require simultaneous temporal landmarks or a specified minimum spatial distance.

## Saved arrays

`temporal_features_with_observed.npz` has 109 named arrays. Element 0 is the observed estimate, and elements 1-20000 are biological-bootstrap draws. Shared draw positions preserve the co-measurement dependence used for paired contrasts.

`temporal_contrast_draws.npz` has 226 named arrays keyed by `contrast_id`, each containing 20,000 B-minus-A bootstrap differences without an observed element. Several Landmark inventory rows use the same contrast array because they evaluate the temporal relationship at different spatial stages.
