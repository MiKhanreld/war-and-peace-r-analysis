library(xml2)
library(tidyverse)
library(igraph)

top_n <- 25

img_width <- 2200
img_height <- 1600
img_res <- 200

seed_value <- 123

node_label_cex <- 0.9
legend_cex <- 0.95

edge_mult <- 2.2
alpha_group <- 0.18

color_palette <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
)

group_names_top <- c(
  "Окружение Болконских",
  "Петербургский свет",
  "Военные",
  "Семья Ростовых",
  "Французская сторона"
)

group_names_vol1 <- c(
  "Болконские",
  "Ростовы",
  "Русские офицеры",
  "Петербургский салон",
  "Военное командование"
)

doc <- read_xml("War_and_Peace.xml")
ns <- xml_ns_rename(xml_ns(doc), d1 = "tei")

clean_ids <- function(x) {
  x |>
    str_replace_all("#", "") |>
    str_replace_all("[^[:alnum:]_ ;]", "_") |>
    str_replace_all("_+", "_") |>
    str_trim()
}

replicas <- xml_find_all(doc, ".//tei:text//tei:div[@type='volume']//tei:said", ns) |>
  map_dfr(function(x) {

    volume_node <- xml_find_first(x, "ancestor::tei:div[@type='volume']", ns)

    tibble(
      volume = xml_attr(volume_node, "n"),
      source = xml_attr(x, "who"),
      target = xml_attr(x, "corresp")
    )

  })

replicas$source <- clean_ids(replicas$source)
replicas$target <- clean_ids(replicas$target)

edges_tbl <- replicas |>
  filter(!is.na(source), !is.na(target)) |>
  separate_rows(target, sep = ";") |>
  mutate(target = str_trim(target)) |>
  filter(source != "", target != "", source != target) |>
  count(volume, source, target, name = "Weight")

write.csv(edges_tbl, "speech_data.csv", row.names = FALSE)

g <- graph_from_data_frame(
  edges_tbl |> select(source, target, Weight),
  directed = TRUE
)

g_u <- as_undirected(
  g,
  mode = "collapse",
  edge.attr.comb = list(Weight = "sum")
)

g_u <- simplify(
  g_u,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = list(Weight = "sum")
)

cat("Сводка по графу\n")
cat("Вершины:", vcount(g), "\n")
cat("Ребра:", ecount(g), "\n")
cat("Плотность:", edge_density(g), "\n")
cat("Компоненты:", components(g)$no, "\n")
cat("Транзитивность:", transitivity(g_u), "\n")

largest_comp <- largest_component(g)

cat("Размер главной компоненты:", vcount(largest_comp), "\n")
cat("Диаметр:", diameter(largest_comp), "\n")

metrics_tbl <- tibble(
  character = V(g)$name,
  degree = degree(g, mode = "all"),
  weighted_degree = strength(g, mode = "all", weights = E(g)$Weight),
  betweenness = betweenness(g_u, normalized = TRUE)[match(V(g)$name, V(g_u)$name)]
) |>
  arrange(desc(weighted_degree))

closeness_tbl <- tibble(
  character = V(largest_comp)$name,
  closeness = closeness(largest_comp, mode = "all", normalized = TRUE)
) |>
  arrange(desc(closeness))

write.csv(metrics_tbl, "top_metrics.csv", row.names = FALSE)
write.csv(closeness_tbl, "closeness_largest_component.csv", row.names = FALSE)

cut_points <- V(g_u)$name[articulation_points(g_u)]

write.csv(
  tibble(character = cut_points),
  "articulation_points.csv",
  row.names = FALSE
)

cat("Размер наибольшей клики:\n")
print(clique_num(g_u))

k_cores <- coreness(g_u)

write.csv(
  tibble(character = V(g_u)$name, core = k_cores),
  "kcore_values.csv",
  row.names = FALSE
)

walk_all <- cluster_walktrap(g_u, weights = E(g_u)$Weight)

groups_all <- membership(walk_all)
mod_all <- modularity(walk_all)

write.csv(
  tibble(character = V(g_u)$name, group = groups_all),
  "groups_all.csv",
  row.names = FALSE
)

cat("Модулярность общего графа:\n")
print(mod_all)

V(g_u)$wdegree <- strength(g_u, weights = E(g_u)$Weight)

top_names <- V(g_u)$name[
  order(V(g_u)$wdegree, decreasing = TRUE)[1:top_n]
]

sub_g <- induced_subgraph(g_u, top_names)

walk_sub <- cluster_walktrap(sub_g, weights = E(sub_g)$Weight)
groups_sub <- membership(walk_sub)

node_colors_sub <- color_palette[
  (groups_sub - 1) %% length(color_palette) + 1
]

group_list_sub <- split(V(sub_g), groups_sub)

png(
  "network_top25.png",
  width = img_width,
  height = img_height,
  res = img_res
)

set.seed(seed_value)

plot(
  walk_sub,
  sub_g,
  layout = layout_with_fr(
    sub_g,
    niter = 5000,
    area = vcount(sub_g)^2.8,
    repulserad = vcount(sub_g)^3.5
  ),
  col = node_colors_sub,
  vertex.frame.color = "black",
  vertex.label = V(sub_g)$name,
  vertex.label.cex = node_label_cex,
  vertex.label.color = "navy",
  vertex.size = 9 + V(sub_g)$wdegree / max(V(sub_g)$wdegree) * 20,
  edge.color = adjustcolor("grey45", alpha.f = 0.7),
  edge.width = E(sub_g)$Weight / max(E(sub_g)$Weight) * edge_mult,
  mark.groups = group_list_sub,
  mark.col = adjustcolor(
    color_palette[unique(groups_sub)],
    alpha.f = alpha_group
  ),
  mark.border = color_palette[unique(groups_sub)],
  main = "Топ-25 персонажей по группам"
)

legend(
  "topright",
  legend = group_names_top[
    seq_along(sort(unique(groups_sub)))
  ],
  fill = color_palette[sort(unique(groups_sub))],
  border = color_palette[sort(unique(groups_sub))],
  bty = "n",
  cex = legend_cex
)

dev.off()

edges_vol1 <- edges_tbl |> filter(volume == "1")

g1 <- graph_from_data_frame(
  edges_vol1 |> select(source, target, Weight),
  directed = TRUE
)

g1_u <- as_undirected(
  g1,
  mode = "collapse",
  edge.attr.comb = list(Weight = "sum")
)

g1_u <- simplify(
  g1_u,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = list(Weight = "sum")
)

V(g1_u)$wdegree <- strength(g1_u, weights = E(g1_u)$Weight)

top_names_1 <- V(g1_u)$name[
  order(V(g1_u)$wdegree, decreasing = TRUE)[1:top_n]
]

sub_g1 <- induced_subgraph(g1_u, top_names_1)

walk_1 <- cluster_walktrap(sub_g1, weights = E(sub_g1)$Weight)

groups_1 <- membership(walk_1)

node_colors_1 <- color_palette[
  (groups_1 - 1) %% length(color_palette) + 1
]

group_list_1 <- split(V(sub_g1), groups_1)

write.csv(
  tibble(
    character = V(sub_g1)$name,
    weighted_degree = V(sub_g1)$wdegree,
    group = groups_1
  ),
  "groups_volume1.csv",
  row.names = FALSE
)

cat("Модулярность тома 1:\n")
print(modularity(walk_1))

png(
  "groups_volume1.png",
  width = img_width,
  height = img_height,
  res = img_res
)

set.seed(seed_value)

plot(
  walk_1,
  sub_g1,
  layout = layout_with_fr(
    sub_g1,
    niter = 5000,
    area = vcount(sub_g1)^2.8,
    repulserad = vcount(sub_g1)^3.5
  ),
  col = node_colors_1,
  vertex.frame.color = "black",
  vertex.label = V(sub_g1)$name,
  vertex.label.cex = node_label_cex,
  vertex.label.color = "navy",
  vertex.size = 9 + V(sub_g1)$wdegree / max(V(sub_g1)$wdegree) * 20,
  edge.color = adjustcolor("grey45", alpha.f = 0.7),
  edge.width = E(sub_g1)$Weight / max(E(sub_g1)$Weight) * edge_mult,
  mark.groups = group_list_1,
  mark.col = adjustcolor(
    color_palette[unique(groups_1)],
    alpha.f = alpha_group
  ),
  mark.border = color_palette[unique(groups_1)],
  main = "Группы в томе 1"
)

legend(
  "topright",
  legend = group_names_vol1[
    seq_along(sort(unique(groups_1)))
  ],
  fill = color_palette[sort(unique(groups_1))],
  border = color_palette[sort(unique(groups_1))],
  bty = "n",
  cex = legend_cex
)

dev.off()