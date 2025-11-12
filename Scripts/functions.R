

# Fonctions pour les définitions des voisinages entre AU ####

library(tidyverse)
library(sf)
library(sp)
library(spdep)
# remotes::install_github("adeverse/adespatial")  # voir : https://github.com/adeverse/adespatial
## tuto : https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html
# install.packages("adespatial")
library(adespatial)


# Fonction pour calculer le voisinage des polygones 


calculate_polygon_neighbors <- function(polygons, dist, symmetric) {
  nb <- st_is_within_distance(polygons, polygons, dist)
  
  # Enlever l'auto-référence et gérer les unités sans voisins
  nb <- lapply(seq_along(nb), function(i) {
    neighbors <- setdiff(nb[[i]], i)
    if (length(neighbors) == 0) {
      return(0L)  # Attribuer 0 aux unités sans voisins
    } else {
      return(as.integer(neighbors))
    }
  })
  
  # Assurer que la liste de voisins a la bonne classe et les bons attributs
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_along(nb))
  attr(nb, "call") <- match.call()
  
  if (symmetric) {
    nb <- try(spdep::make.sym.nb(nb), silent = TRUE)
    if (inherits(nb, "try-error")) {
      return(list(nb = nb, nb_obj = NULL, error = "Erreur lors de la création de la liste symétrique"))
    }
  }
  
  nb_obj <- try(nb2listw(nb, style = "B", zero.policy = TRUE), silent = TRUE)
  if (inherits(nb_obj, "try-error")) {
    return(list(nb = nb, nb_obj = NULL, error = "Erreur lors de la création de l'objet listw"))
  }
  
  list(nb = nb, nb_obj = nb_obj, error = NULL)
}


# Fonction pour calculer le voisinage des centroïdes


calculate_centroid_neighbors <- function(centroids, method, param, symmetric) {
  coords <- st_coordinates(centroids)
  if (method == "knn") {
    # Ajuster k si nécessaire
    n_points <- nrow(coords)
    adjusted_k <- min(param, n_points - 1)  # k ne peut pas être plus grand que n-1
    
    if (adjusted_k < param) {
      warning(sprintf("k ajusté de %d à %d en raison du nombre limité de points", param, adjusted_k))
    }
    
    if (adjusted_k > 0) {
      knn <- knearneigh(coords, k = adjusted_k)
      nb <- knn2nb(knn, sym = symmetric)
    } else {
      # Si adjusted_k est 0 (cas avec 1 seul point)
      nb <- list(0L)
      class(nb) <- "nb"
    }
  } else if (method == "dist") {
    nb <- dnearneigh(coords, d1 = 0, d2 = param)
    if (symmetric) {
      nb <- make.sym.nb(nb)
    }
  } else if (method == "gabriel") {
    nb <- graph2nb(gabrielneigh(coords), sym = symmetric)
  }
  
  # Vérifier si la liste de voisins est vide
  if (all(sapply(nb, length) == 0)) {
    return(list(nb = nb, nb_obj = NULL, error = "Aucun voisin trouvé"))
  }
  
  nb_obj <- try(nb2listw(nb, style = "B", zero.policy = TRUE), silent = TRUE)
  if (inherits(nb_obj, "try-error")) {
    return(list(nb = nb, nb_obj = NULL, error = "Erreur lors de la création de l'objet listw"))
  }
  
  list(nb = nb, nb_obj = nb_obj, error = NULL)
}

# Fonction pour créer toutes les configurations de voisinage


create_neighborhood_configs <- function(polygons, centroids, k_values, dist_values, dist_values_points, group_column = NULL) {
  # Configurations pour les polygones (uniquement distance)
  polygon_configs <- expand.grid(
    type = "polygon",
    method = "dist",
    param = dist_values,
    symmetric = FALSE
  )
  
  # Configurations pour les centroïdes (knn et distance)
  centroid_configs <- rbind(
    expand.grid(
      type = "centroid",
      method = "knn",
      param = k_values,
      symmetric = c(FALSE, TRUE)
    ),
    expand.grid(
      type = "centroid",
      method = "dist",
      param = dist_values_points,
      symmetric = FALSE
    ),
    data.frame(
      type = "centroid",
      method = "gabriel",
      param = NA,
      symmetric = TRUE
    )
  )
  
  # Combiner toutes les configurations
  configs <- rbind(polygon_configs, centroid_configs)
  
  # Fonction pour calculer les voisinages pour un groupe donné
  calculate_group_neighbors <- function(group_polygons, group_centroids, config) {
    result <- if (config$type == "polygon") {
      calculate_polygon_neighbors(group_polygons, config$param, config$symmetric)
    } else {
      calculate_centroid_neighbors(group_centroids, config$method, 
                                   ifelse(is.na(config$param), NULL, config$param), 
                                   config$symmetric)
    }
    return(result)
  }
  
  # Traitement avec groupes régionaux
  neighborhood_list <- map(1:nrow(configs), function(i) {
    config <- configs[i,]
    
    if (!is.null(group_column)) {
      # Obtenir les groupes uniques
      groups <- unique(st_drop_geometry(polygons)[[group_column]])
      
      # Créer une liste de voisinage de la taille totale des données
      n_total <- nrow(polygons)
      combined_nb <- vector("list", n_total)
      
      # Traiter chaque groupe séparément
      for (group in groups) {
        # Filtrer les données pour le groupe actuel
        group_mask <- st_drop_geometry(polygons)[[group_column]] == group
        group_polygons <- polygons[group_mask, ]
        group_centroids <- centroids[group_mask, ]
        
        # Obtenir les indices originaux pour ce groupe
        original_indices <- which(group_mask)
        
        # Calculer les voisinages pour ce groupe
        group_result <- if (config$type == "polygon") {
          calculate_polygon_neighbors(group_polygons, config$param, config$symmetric)
        } else {
          calculate_centroid_neighbors(group_centroids, config$method, 
                                       ifelse(is.na(config$param), NULL, config$param), 
                                       config$symmetric)
        }
        
        if (!is.null(group_result$error)) {
          next
        }
        
        # Convertir les indices régionaux en indices globaux
        for (j in seq_along(group_result$nb)) {
          if (length(group_result$nb[[j]]) > 0 && group_result$nb[[j]][1] != 0) {
            local_indices <- group_result$nb[[j]]
            global_indices <- original_indices[local_indices]
            combined_nb[[original_indices[j]]] <- global_indices
          } else {
            combined_nb[[original_indices[j]]] <- 0L
          }
        }
      }
      
      # S'assurer que les positions sans voisins ont 0L
      combined_nb[sapply(combined_nb, is.null)] <- 0L
      
      # Ajouter les attributs nécessaires
      class(combined_nb) <- "nb"
      attr(combined_nb, "region.id") <- as.character(1:length(combined_nb))
      attr(combined_nb, "call") <- match.call()
      
      # Créer l'objet listw
      combined_nb_obj <- try(nb2listw(combined_nb, style = "B", zero.policy = TRUE), silent = TRUE)
      if (inherits(combined_nb_obj, "try-error")) {
        return(list(config = config, error = "Erreur lors de la création de l'objet listw"))
      }
      
      result <- list(
        nb = combined_nb,
        nb_obj = combined_nb_obj
      )
      
    } else {
      # Calcul sans groupes (comme avant)
      result <- if (config$type == "polygon") {
        calculate_polygon_neighbors(polygons, config$param, config$symmetric)
      } else {
        calculate_centroid_neighbors(centroids, config$method, 
                                     ifelse(is.na(config$param), NULL, config$param), 
                                     config$symmetric)
      }
      
      if (!is.null(result$error)) {
        return(list(config = config, error = result$error))
      }
    }
    
    # Calculer les statistiques sur l'objet final
    list(
      config = config,
      nb = result$nb,
      nb_obj = result$nb_obj,
      n_components = n.comp.nb(result$nb)$nc,
      avg_neighbors = mean(card(result$nb)),
      med_neighbors = median(card(result$nb)),
      min_neighbors = min(card(result$nb)),
      max_neighbors = max(card(result$nb)),
      q1_neighbors = quantile(card(result$nb), prob = 0.25),
      q3_neighbors = quantile(card(result$nb), prob = 0.75)
    )
  })
  
  # Filtrer les configurations qui ont échoué
  valid_configs <- Filter(function(x) is.null(x$error), neighborhood_list)
  
  # Nommer les éléments de la liste
  named_configs <- lapply(valid_configs, function(config) {
    name <- paste(config$config$type, 
                  config$config$method, 
                  ifelse(is.na(config$config$param), "NA", config$config$param), 
                  config$config$symmetric, 
                  sep = "_")
    setNames(list(config), name)
  })
  
  # Combiner la liste nommée
  result <- do.call(c, named_configs)
  
  return(result)
}



# Fonction de visualisation 
create_simple_plot <- function(config, geometry) {
  title_text <- paste(
    "Type:", config$config$type,
    "| Méthode:", config$config$method,
    "| Paramètre:", config$config$param,
    "\nSymétrique:", ifelse(config$config$symmetric, "Oui", "Non"),
    "| Nb composantes:", config$n_components,
    "| Nb moyen voisins:", round(config$avg_neighbors, 2)
    
  )
  
  plot(st_geometry(geometry), pch = 20, col = "black", main = title_text)
  plot(config$nb, st_coordinates(st_centroid(geometry)), add = TRUE, col = "red", lwd = 0.5)
  
}



# Fonction pour le calcul des MORAN globaux ####


# Calcul des tests de Moran NP pour chaque couple variable-liste de voisins
calculate_moran_tests <- function(data, neighborhood_configs) {
  variables <- select(data, where(is.numeric))
  
  results <- map(neighborhood_configs, function(config) {
    if (is.null(config$nb_obj)) {
      return(NULL)
    }
    
    var_results <- map(names(variables), function(var_name) {
      tryCatch({
        moranNP <- moranNP.randtest(variables[[var_name]], config$nb_obj, nrepet = 999, alter = "two-sided")
        
        data.frame(
          variable = var_name,
          type = c("I+", "I-"),
          observed = moranNP$obs,
          pvalue = moranNP$pvalue,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        warning(paste("Error in Moran NP test for variable:", var_name))
        NULL
      })
    })
    
    var_results <- bind_rows(var_results)
    
    if (nrow(var_results) > 0) {
      var_results %>%
        mutate(
          config_type = config$config$type,
          method = config$config$method,
          param = config$config$param,
          symmetric = config$config$symmetric,
          n_components = config$n_components
        )
    } else {
      NULL
    }
  })
  
  bind_rows(results)
}


# FONCTION MORAN LOCAUX ####


# exploration systématique du nombre de moran significatifs
analyze_local_moran <- function(data, var_names, listw_configs, seed = 123, seuil) {
  # Pour la reproductibilité
  set.seed(seed)
  
  # Création d'une liste pour stocker les résultats résumés
  results_summary <- list()
  # Création d'une liste pour stocker les résultats détaillés
  detailed_results <- list()
  
  # Pour chaque variable
  for(var_name in var_names) {
    # Vérification que la variable existe dans le dataframe
    if(!var_name %in% names(data)) {
      warning(paste("Variable", var_name, "non trouvée dans le jeu de données. Ignorée."))
      next
    }
    
    # Pour chaque configuration de voisinage
    for(listw_name in names(listw_configs)) {
      # Calcul du Moran local
      lm_result <- localmoran_perm(
        data[[var_name]], 
        listw = listw_configs[[listw_name]]$nb_obj, 
        nsim = 999
      )
      
      # Création du dataframe avec les résultats
      df_result <- as.data.frame(bind_cols(lm_result, attr(lm_result, "quadr"))) %>%
        mutate(
          p.value.holm = p.adjust(`Pr(z != E(Ii))`, method = "holm"),
          p.value.fdr = p.adjust(`Pr(z != E(Ii))`, method = "fdr"),
          p.value.bonferroni = p.adjust(`Pr(z != E(Ii))`, method = "bonferroni")
        )
      
      # Stockage du résultat détaillé
      result_key <- paste(var_name, listw_name, sep = "_")
      detailed_results[[result_key]] <- df_result
      
      # Comptage des significatifs pour chaque méthode
      sig_counts <- list(
        holm = sum(df_result$p.value.holm < seuil, na.rm = TRUE),
        fdr = sum(df_result$p.value.fdr < seuil, na.rm = TRUE),
        bonferroni = sum(df_result$p.value.bonferroni < seuil, na.rm = TRUE)
      )
      
      # Comptage des modalités pour les cas significatifs (FDR)
      sig_cases <- df_result[df_result$p.value.fdr < seuil, ]
      
      # Création d'un vecteur avec tous les comptages possibles initialisés à 0
      quad_counts <- c(
        "High-High" = 0,
        "Low-Low" = 0,
        "Low-High" = 0,
        "High-Low" = 0
      )
      
      # Mise à jour avec les valeurs réelles seulement s'il y a des cas significatifs
      if(nrow(sig_cases) > 0) {
        quadrant_counts <- table(sig_cases$mean)
        quad_counts[names(quadrant_counts)] <- quadrant_counts
      }
      
      # Stockage des résultats résumés
      results_summary[[result_key]] <- c(
        variable = var_name,
        voisinage_km = parse_number(listw_name)/1000,
        sig_holm = sig_counts$holm,
        sig_fdr = sig_counts$fdr,
        sig_bonferroni = sig_counts$bonferroni,
        quad_counts
      )
    }
  }
  
  # Conversion du résumé en dataframe
  results_summary_df <- bind_rows(lapply(results_summary, function(x) as.data.frame(t(x)))) %>%
    mutate(across(c(sig_holm, sig_fdr, sig_bonferroni, 
                    `High-High`, `Low-Low`, `Low-High`, `High-Low`), 
                  ~as.numeric(replace(., is.na(.), 0))))
  
  # Retourne une liste avec les deux types de résultats
  return(list(
    summary = results_summary_df,
    detailed = detailed_results
  ))
}
