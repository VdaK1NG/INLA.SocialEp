#' Make a fully connected graph for INLA
#'
#' This function converts a graph with islands into a fully connected graph
#' by searching for the nearest area between the different nodes and connecting
#' them, with the possibility of setting a threshold for connecting several areas
#' if they are within the limits.
#'
#' @param sp_obj Spatial object in sp form
#' @param graph Graph object from nb2INLA function
#' @param snap Snap argument for the function poly2nb
#' @param distance_by Choices are perimeter and centroid
#' @return List including three different objects: the adjacency object, the edited graph and a dataframe with the changes applied
#' @export

inla.full.graph <- function(sp_obj, graph, thresh = 0.10, snap = 5e-07, distance_by = c("perimeter", "centroid")) {

  # Detect distance metric
  distance_by <- match.arg(distance_by)
  sp_obj_sf <- st_as_sf(sp_obj)

  # Basic checks
  if (is.null(graph$cc) || is.null(graph$cc$nodes)) stop("Graph does not have the expected structure.")

  # Turn thresh into the correct value
  if (!is.numeric(thresh) || length(thresh) != 1 || thresh < 0) stop("Thresh must be a number >= 0 (ex. 0.1 or 10).")
  thresh_prop <- if (thresh > 1) thresh / 100 else thresh

  n_nodes <- graph$cc$nodes
  if (length(n_nodes) == 1) {
    cat("-> Graph only has 1 node, no need to fix it\n")
    return(list(adj = poly2nb(sp_obj, snap = snap), graph = graph, connected_pairs = NULL))
  } else {
    cat(paste0("-> Graph has ", length(n_nodes), " components; connecting with thresh = ", thresh_prop, "\n"))
  }

  # Keep geometry only to avoid future warnings
  geom_only <- st_geometry(sp_obj_sf)

  # Choose distance metrix
  if (distance_by == "perimeter") {
    boundary_geom <- st_boundary(geom_only)
    is_empty <- st_is_empty(boundary_geom)
    if (any(is_empty)) {
      warning("Some geomteries are empy; using centroid as backup for this areas.")
      cent_geom <- st_centroid(geom_only)
      boundary_geom[is_empty] <- cent_geom[is_empty]
    }
    mtx_raw <- st_distance(boundary_geom, boundary_geom)
  } else {
    cent_geom <- st_centroid(geom_only)
    mtx_raw <- st_distance(cent_geom, cent_geom)
  }

  # pick stable ids (prefer rownames(sp_obj), else 1..n)
  rid <- rownames(sp_obj_sf)
  if (is.null(rid)) rid <- as.character(seq_len(nrow(sp_obj_sf)))

  # numeric matrix (sin unidades)
  mtx_distance <- matrix(as.numeric(mtx_raw),
                         nrow = nrow(mtx_raw), ncol = ncol(mtx_raw),
                         dimnames = list(rid, rid))

  # adyacencia base (lista nb) y copia segura para editar
  adj_orig <- poly2nb(sp_obj, snap = snap)

  # copiar adj_orig en adj_ed asegurando tipo integer y sin modificar elementos existentes
  adj_ed <- adj_orig

  # tabla para auditar conexiones añadidas
  connected <- list()

  # componente base (empezamos por el primer componente)
  areas_n1 <- n_nodes[[1]]

  for (i in seq(2, length(n_nodes))) {
    areas_n2 <- n_nodes[[i]]

    if (is.numeric(areas_n1)) areas_n1 <- rid[areas_n1]
    if (is.numeric(areas_n2)) areas_n2 <- rid[areas_n2]

    # areas_n1 / areas_n2 might be indices; convert to names for safe subsetting
    areas_n1 <- as.character(areas_n1)
    areas_n2 <- as.character(areas_n2)

    if (length(areas_n2) == 0) next

    # submatriz (filas = areas_n1, cols = areas_n2)
    subm <- mtx_distance[areas_n1, areas_n2, drop = FALSE]

    # si todas NA -> saltar y expandir base
    if (all(is.na(subm))) {
      warning(sprintf("---> All distances NA between base component and component %d. Skipping.", i))
      areas_n1 <- unique(c(areas_n1, areas_n2))
      next
    }

    # distancia mínima entre conjuntos
    dist_min <- min(subm, na.rm = TRUE)
    # límite según thresh relativo
    limit <- dist_min * (1 + thresh_prop)

    # seleccionar todos los pares cuya distancia <= limit
    pairs <- which(!is.na(subm) & subm <= limit, arr.ind = TRUE)


    # IMPORTANT: adj_ed is indexed by INTEGER area ids (1..N), so we must map row/col to those integers.
    # We have subm rownames/colnames as rid labels, so map them back to indices:
    rid_to_idx <- setNames(seq_along(rid), rid)

    for (r in seq_len(nrow(pairs))) {
      row_idx <- pairs[r, "row"]   # position inside areas_n1
      col_idx <- pairs[r, "col"]   # position inside areas_n2

      # these are the rid labels used in subm dimnames
      node1_name <- areas_n1[row_idx]
      node2_name <- areas_n2[col_idx]

      # convert back to numeric indices for nb object
      node1 <- as.integer(rid_to_idx[[node1_name]])
      node2 <- as.integer(rid_to_idx[[node2_name]])

      dist_val <- subm[row_idx, col_idx]

      if (is.na(node1) || is.na(node2)) {
        warning(sprintf("---> Could not map ids to indices: %s - %s. Skipping.", node1_name, node2_name))
        next
      }

      adj_ed[[node1]] <- c(adj_ed[[node1]], node2)
      adj_ed[[node2]] <- c(adj_ed[[node2]], node1)

      connected[[length(connected) + 1]] <- data.frame(
        node1 = node1,
        node2 = node2,
        node1_id = node1_name,
        node2_id = node2_name,
        distance = as.numeric(dist_val),
        dist_min = as.numeric(dist_min),
        limit = as.numeric(limit),
        component_joined = i,
        stringsAsFactors = FALSE
      )
      cat(sprintf("---> (added) Connected %d(%s) <-> %d(%s) (d=%g, dist_min=%g, limit=%g)\n",
                  node1, node1_name, node2, node2_name, dist_val, dist_min, limit))
    }

    # expandir la base con el componente ya unido (sin modificar adj_ed de otros nodos)
    areas_n1 <- unique(c(areas_n1, areas_n2))
  }

  # convertir registros a data.frame (si hay)
  connected_pairs <- if (length(connected) == 0) {
    data.frame(node1 = integer(0), node2 = integer(0), distance = numeric(0),
               dist_min = numeric(0), limit = numeric(0), component_joined = integer(0))
  } else {
    do.call(rbind, connected)
  }

  # Create INLA graph and remove temporal file
  grafo_new <- NULL
  tryCatch({
    nb2INLA("MapGraph", adj_ed)
    grafo_new <- inla.read.graph("MapGraph")
    files_to_remove <- list.files(pattern = paste0("^", "MapGraph"))
    if (length(files_to_remove) > 0) file.remove(files_to_remove)
    cat("-> INLA graph written successfully.\n")
  }, error = function(e) {
    stop("Error while writing INLA graph: ", conditionMessage(e))
  })

  return(list(adj = adj_ed, graph = grafo_new, connected_pairs = connected_pairs))
}
