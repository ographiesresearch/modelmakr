library(sf)
library(elevatr)
library(raster)
library(dplyr)
library(tigris)
library(units)
library(lwgeom)
library(stringr)
library(purrr)
library(osmdata)

options(tigris_use_cache = TRUE)

elevate <- function(sdf) {
  crs <- sf::st_crs(sdf)
  sdf %>%
    dplyr::mutate(
      wkt = lwgeom::st_asewkt(geometry),
      wkt = stringr::str_replace_all(wkt, '(?:[0-9]),', paste0(" ", z, ",")),
      wkt = stringr::str_replace_all(wkt, '(?:[0-9])\\)', paste0(" ", z, ")"))
    ) %>% 
    sf::st_drop_geometry() %>%
    sf::st_as_sf(wkt = "wkt", crs = crs)
}

model_scale <- function(sdf, scale_elements, thickness, z_scale, contour_interval) {
  crs <- sf::st_crs(sdf)
  sdf_scaled <- (sf::st_geometry(sdf) - c(scale_elements['xadjust'], scale_elements['yadjust'])) * scale_elements['scale_factor']
  
  if (!("sfc" %in% class(sdf))) {
    sf::st_geometry(sdf) <- sdf_scaled
    if ("z" %in% colnames(sdf))
      sdf <- sdf %>% 
        dplyr::mutate(z = (z / (z_scale * contour_interval)) * as.numeric(set_units(thickness, m)))
    print(sdf)
  }
  sdf <- sdf %>%
    sf::st_set_crs(crs)
  sdf
}

scale_elements <- function(sdf, model_size = NULL) {
  if (!is.null(model_size)) {
    bbox <- sf::st_bbox(sdf)
    xsize <- abs(bbox$xmax - bbox$xmin)
    ysize <- abs(bbox$ymax - bbox$ymin)
    xcenter <- as.numeric(bbox$xmax - (xsize / 2))
    ycenter <- as.numeric(bbox$ymax - (ysize / 2))
    
    largestdim <- max(xsize, ysize)
    scale_factor <- as.numeric(set_units(model_size, m)) / largestdim
  } else {
    xcenter <- 0
    ycenter <- 0
    scale_factor <- 1
  }
  c(scale_factor = scale_factor, xadjust = xcenter, yadjust = ycenter)
}

write_to_formats <- function(sdf, what, out_dir, exts, scale_elements, thickness = NULL, z_scale, contour_interval) {
  if ("geojson" %in% exts) {
    sf::st_write(
      sf::st_transform(sdf, 4326), 
      base::file.path(out_dir, base::paste(what, "geojson", sep = ".")),
      driver = "geojson", 
      append = FALSE, 
      delete_dsn = TRUE
    )
  }
  
  if ("dxf" %in% exts) {
    sdf <- sdf %>%
      model_scale(scale_elements = scale_elements, thickness = thickness, z_scale, contour_interval)
    
    if ("z" %in% colnames(sdf)) {
      sdf <- elevate(sdf)
    }
    
    sdf %>%
      sf::st_geometry() %>%
      sf::st_write(
        base::file.path(out_dir, base::paste(what, "dxf", sep = ".")),
        driver = "dxf", 
        append = FALSE, 
        delete_dsn = TRUE
      )
  }
  # Return original SF dataframe.
  sdf
}

get_crs <- function(center, util_dir, base_crs) {
  utm_zone <- center %>%
    sf::st_join(
      sf::st_read(
        base::file.path(util_dir, base::paste("utm_zones", "geojson", sep = "."))
      ) %>%
        st_transform(base_crs),
      sf::st_intersects,
      largest = TRUE
    ) %>%
    dplyr::pull(zone_num)
  
  center_coords <- center %>%
    sf::st_transform(4326) %>%
    sf::st_coordinates()
  
  # UTM convention is that unmarked zones are
  # N, S zones are marked S (because Eurocentrism).
  # Here, we test if latitude > 0 to determine
  # whether it is N or S.
  if(center_coords[1,2] < 0){
    utm_zone <- base::paste0(utm_zone, "S")
  }
  
  # Set up projection using custom proj4string.
  sf::st_crs(
    base::paste0(
      "+proj=utm +zone=", 
      utm_zone, 
      " +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
    )
  )
}

get_basemap <- function(
    area = 'data/Neighborhoods.shp',
    buffer = 1.25,
    year = 2020,
    exts = c("geojson", "dxf"),
    out_dir = "results",
    shape = "c", # TODO: or "s"
    z_scale = 4,
    contour_interval = 3,
    util_dir = "data",
    base_crs = 3857,
    model_size = as_units(6, "ft"),
    thickness = NULL
){
  if (!file.exists(out_dir)) {
    dir.create(file.path(getwd(), out_dir))
  }
  if (is.null(thickness)) {
    thickness <- contour_interval
  }
  
  if (!base::is.null(area)) {
    shape <- sf::st_read(area) %>%
      sf::st_union() %>%
      sf::st_transform(base_crs)
    
    center <- shape %>%
      sf::st_point_on_surface() %>%
      sf::st_as_sf() 
    
    dist <- center %>%
      sf::st_distance(sf::st_cast(shape, "POINT")) %>%
      base::max()
    
    dl_radius <- center %>%
      sf::st_buffer(dist * (buffer ^ 2))
    
    study_radius <- center %>%
      sf::st_buffer(dist * buffer)
    
  } else if (!base::is.null(coords)) {
    
  } else {
    base::print("Expected either an area or a coordinate.")
    break
  }
  
  # Automatically get UTM zone projection.
  crs <- get_crs(center, util_dir, base_crs)
  
  # transform study area radii.
  dl_radius <- dl_radius %>%
    sf::st_transform(crs)
  
  study_radius <- study_radius %>%
    sf::st_transform(crs)
  
  se <- scale_elements(study_radius, model_size = model_size)
  
  study_radius %>%
    write_to_formats(
      "study_radius", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  states <- states(year = year) %>%
    sf::st_transform(crs) %>%
    sf::st_filter(
      study_radius,
      .predicate = sf::st_intersects
    ) %>%
    dplyr::pull(STATEFP)
  
  counties <- tigris::counties(states, year = year) %>%
    sf::st_transform(crs) %>%
    sf::st_filter(
      study_radius,
      .predicate = sf::st_intersects
    ) %>%
    dplyr::pull(COUNTYFP)
  
  water <- tigris::area_water(
    state = states,
    county = counties,
    year = year
  ) %>%
    sf::st_transform(crs) %>%
    sf::st_intersection(dl_radius) %>%
    dplyr::mutate(
      area = sf::st_area(geometry)
    ) %>%
    dplyr::filter(
      area > units::set_units(25000, m^2)
    ) %>%
    sf::st_union() %>%
    st_as_sf()
  
  dem <- elevatr::get_elev_raster(
      locations = st_transform(dl_radius, 4326),
      z = 14,
      clip = "locations"
    ) %>%
    projectRaster(crs = crs$input) %>%
    reclassify(cbind(-Inf, 0, 0), right = FALSE)
  
  
  dem %>%
    raster::mask(study_radius) %>%
    raster::writeRaster(
      base::file.path(out_dir, base::paste("dem", "tif", sep = ".")),
      overwrite = TRUE
    )
  
  contours <- dem %>%
    raster::rasterToContour(
      levels = base::seq(
        from = base::floor(raster::minValue(.)),
        to = base::ceiling(raster::maxValue(.)),
        by = contour_interval
      )
    ) %>%
    sf::st_as_sf(base::seq_len(base::length())) %>%
    dplyr::rename(
      z = level
    ) %>%
    sf::st_cast("LINESTRING") %>%
    dplyr::mutate(
      z = base::as.integer(z) * z_scale,
      length = sf::st_length(geometry)
    ) %>%
    dplyr::filter(
      length > units::set_units(250, m)
    )
  
  contour_polys <- study_radius %>%
    lwgeom::st_split(contours) %>%
    sf::st_collection_extract("POLYGON") %>%
    sf::st_set_geometry("geometry") %>%
    tibble::rowid_to_column("id")
  
  poly_to_contour <- contour_polys %>%
    sf::st_join(contours) %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(
      z_min = base::as.integer(base::min(z)),
      z_max = base::as.integer(base::max(z)),
      z = base::min(z_min, z_max)
    )
  
  contour_joined <- contour_polys %>%
    dplyr::left_join(poly_to_contour, by = c("id" = "id")) %>%
    dplyr::group_by(z) %>%
    dplyr::summarize(
      geometry = sf::st_cast(sf::st_union(geometry), "MULTIPOLYGON")
    ) %>%
    arrange(desc(z))
  
  # Export enclosed contours.
  contour_joined %>%
    dplyr::mutate(
      geometry = sf::st_as_sfc(purrr::accumulate(geometry, sf::st_union), crs = crs)
    ) %>%
    sf::st_cast("MULTILINESTRING") %>%
    write_to_formats(
      "contours", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  # Export above contours for alignment.
  contour_joined  %>%
    dplyr::mutate(
      z = lead(z)
    ) %>%
    sf::st_cast("MULTILINESTRING") %>%
    tidyr::drop_na(z) %>%
    st_set_geometry("geometry") %>%
    write_to_formats(
      "contour_guides", 
      out_dir, exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  tigris::roads(
    state = states,
    county = counties,
    year = year
  ) %>%
    sf::st_transform(crs) %>%
    sf::st_intersection(contour_joined) %>%
    write_to_formats("roads", out_dir, exts, scale_elements = se, thickness = thickness, z_scale = z_scale, contour_interval = contour_interval)
  
  tigris::rails(year = year) %>%
    sf::st_transform(crs) %>%
    sf::st_intersection(contour_joined) %>%
    write_to_formats(
      "rails", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  shape %>%
    sf::st_transform(crs) %>%
    sf::st_cast("MULTILINESTRING") %>%
    sf::st_intersection(contour_joined) %>%
    write_to_formats(
      "shape", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  shape %>%
    sf::st_transform(crs) %>%
    write_to_formats(
      "shape_poly", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  # buildings
  opq(st_bbox(st_transform(dl_radius, 4326))) %>%
    add_osm_feature(key = 'building') %>%
    osmdata_sf() %>%
    extract2("osm_polygons") %>%
    # sf::st_cast("MULTILINESTRING") %>%
    sf::st_transform(crs) %>%
    sf::st_intersection(contour_joined) %>%
    st_set_geometry("geometry") %>%
    write_to_formats(
      "buildings", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  water %>%
    sf::st_intersection(contour_joined) %>%
    sf::st_cast("MULTILINESTRING") %>%
    st_set_geometry("geometry") %>%
    write_to_formats(
      "water", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  
  water %>%
    sf::st_intersection(contour_joined) %>%
    st_set_geometry("geometry") %>%
    write_to_formats(
      "water_poly", 
      out_dir, 
      exts, 
      scale_elements = se, 
      thickness = thickness, 
      z_scale = z_scale, 
      contour_interval = contour_interval
    )
  return(study_radius)
}

trim_other <- function(study_area, geos, name, out_dir = "results") {
  st_read(geos) %>%
    st_transform(st_crs(study_area)) %>%
    st_intersection(st_geometry(study_area)) %>%
    st_transform(4326) %>%
    st_write(
      base::file.path(out_dir, base::paste(name, "geojson", sep = ".")),
      delete_dsn = TRUE)
}

# study_area <- get_basemap(area = "data/neighborhood.geojson", z_scale = 6, model_size = NULL)


# trim_other(
#   study_area = study_area,
#   geos = "/Users/ehuntley/Desktop/dev/eviction_exhibit/data/openspace/OPENSPACE_POLY.shp",
#   name = "open_space"
# )