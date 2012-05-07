
// STL
#include <iostream>
#include <stack>
#include <set>
#include <fstream>

// OGR
#include <gdal/ogrsf_frmts.h>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag PT;
typedef CGAL::Exact_intersections_tag IT;
typedef K::Point_3   Point;
//typedef Triangulation::Point Point;

typedef CGAL::Projection_traits_xy_3<K>  Kxy;
typedef CGAL::Triangulation_vertex_base_2<Kxy> VBxy;
typedef CGAL::Constrained_triangulation_face_base_2<Kxy> FBxy;
typedef CGAL::Triangulation_face_base_with_info_2<void *, Kxy, FBxy> FBWIxy;
typedef CGAL::Triangulation_data_structure_2<VBxy, FBWIxy> TDSxy;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kxy, TDSxy, PT> CDTxy;
typedef CGAL::Constrained_triangulation_plus_2<CDTxy> Triangulationxy;

typedef CGAL::Projection_traits_xz_3<K>  Kxz;
typedef CGAL::Triangulation_vertex_base_2<Kxz> VBxz;
typedef CGAL::Constrained_triangulation_face_base_2<Kxz> FBxz;
typedef CGAL::Triangulation_face_base_with_info_2<void *, Kxz, FBxz> FBWIxz;
typedef CGAL::Triangulation_data_structure_2<VBxz, FBWIxz> TDSxz;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kxz, TDSxz, PT> CDTxz;
typedef CGAL::Constrained_triangulation_plus_2<CDTxz> Triangulationxz;

typedef CGAL::Projection_traits_yz_3<K>  Kyz;
typedef CGAL::Triangulation_vertex_base_2<Kyz> VByz;
typedef CGAL::Constrained_triangulation_face_base_2<Kyz> FByz;
typedef CGAL::Triangulation_face_base_with_info_2<void *, Kyz, FByz> FBWIyz;
typedef CGAL::Triangulation_data_structure_2<VByz, FBWIyz> TDSyz;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kyz, TDSyz, PT> CDTyz;
typedef CGAL::Constrained_triangulation_plus_2<CDTyz> Triangulationyz;


//-- prototypes
bool  triangulateandtag_xy(OGRGeometry* geometry, Triangulationxy &triangulation);
void  tag_xy(Triangulationxy &triangulation, void *interiorHandle, void *exteriorHandle);
bool  triangulateandtag_xz(OGRGeometry* geometry, Triangulationxz &triangulation);
void  tag_xz(Triangulationxz &triangulation, void *interiorHandlez, void *exteriorHandle);
bool  triangulateandtag_yz(OGRGeometry* geometry, Triangulationyz &triangulation);
void  tag_yz(Triangulationyz &triangulation, void *interiorHandle, void *exteriorHandle);

int   get_projection_plane(OGRGeometry* geometry);




int main (int argc, const char * argv[]) {
  
  if (argc < 2 || argc > 3 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    std::cout << "=== prepair Help ===\n" << std::endl;
    std::cout << "Usage:   triface 'POLYGON(...)'" << std::endl;
    std::cout << "OR" << std::endl;
    std::cout << "Usage:   triface -f infile.txt (infile.txt must contain one WKT on the 1st line)" << std::endl;
    return 0;
  }
  
  // Read input
  unsigned int bufferSize = 10000000;
  char *inputWKT = (char *)malloc(bufferSize*sizeof(char *));
  
  for (int argNum = 1; argNum < argc; ++argNum) {
    if (strcmp(argv[argNum], "-f") == 0) {
      
      if (argNum + 1 <= argc - 1 && argv[argNum+1][0] != '-') {
        std::ifstream infile(argv[argNum+1], std::ifstream::in);
        infile.getline(inputWKT, bufferSize);
        ++argNum;
      } else {
        std::cerr << "Error: Missing input file name." << std::endl;
        return 1;
      }
    }
    else 
      strcpy(inputWKT, argv[argNum]);
  }
  
//  std::cout << "Processing: " << inputWKT << std::endl;
  
  OGRGeometry *geometry;
  OGRGeometryFactory::createFromWkt(&inputWKT, NULL, &geometry);
  if (geometry == NULL) {
    std::cout << "Error: WKT is not valid" << std::endl;
    return 1;
  }

  if (geometry->getGeometryType() != wkbPolygon25D) {
    std::cout << "Error: input geometry is not a 3D polygon" << std::endl;
    return 1;
  }
  
  //-- project to proper plane + get flattened geometry
  int proj = get_projection_plane(geometry);
  OGRGeometry *flatgeom = geometry->clone();
  if (proj == 1) {
    OGRPolygon *polygon = (OGRPolygon *)flatgeom;
    for (int curp = 0; curp < polygon->getExteriorRing()->getNumPoints(); ++curp) 
      polygon->getExteriorRing()->setPoint(curp, polygon->getExteriorRing()->getX(curp), polygon->getExteriorRing()->getZ(curp), 0);
  }
  else if (proj == 0) {
    OGRPolygon *polygon = (OGRPolygon *)geometry;
    for (int curp = 0; curp < polygon->getExteriorRing()->getNumPoints(); ++curp)
      polygon->getExteriorRing()->setPoint(curp, polygon->getExteriorRing()->getY(curp), polygon->getExteriorRing()->getZ(curp), 0);
  }
  flatgeom->flattenTo2D();
//  std::cout << "geom: " << geometry->getCoordinateDimension() << std::endl;
//  std::cout << "flatgeom: " << flatgeom->getCoordinateDimension() << std::endl;
  
  //-- check if flattened geometry is valid
  if (flatgeom->IsValid() == FALSE) {
    std::cout << "Error: input polygon is not valid." << std::endl;
    return 1;
  }
  
  if (proj == 2) {
    Triangulationxy triangulation;
    triangulateandtag_xy(geometry, triangulation);
  }
  else if (proj == 1) {
    Triangulationxz triangulation;
    triangulateandtag_xz(geometry, triangulation);
  }
  else { //-- proj == 0
    Triangulationyz triangulation;
    triangulateandtag_yz(geometry, triangulation);    
  }
     
  
  return 0;
//  
//  OGRMultiPolygon* outputPolygons = repair(geometry);
//  
//  if (outputPolygons == NULL) {
//    std::cout << "Impossible to repair the polygon: input points are collinear (no area given)." << std::endl;
//    return 0;
//  }
//  else {
//    char *outputWKT;
//    outputPolygons->exportToWkt(&outputWKT);
//    std::cout << std::endl << "Repaired polygon:" << std::endl << outputWKT << std::endl;
//    return 0;
//  }
}


int get_projection_plane(OGRGeometry* geometry) {
  OGREnvelope3D envelope;
  geometry->getEnvelope(&envelope);
  //-- xy plane
  int proj = 2;
  double projarea = (envelope.MaxX - envelope.MinX) * (envelope.MaxY - envelope.MinY);
  //-- yz plane
  double temp = (envelope.MaxX - envelope.MinX) * (envelope.MaxZ - envelope.MinZ);
  if ( temp > projarea) {
    proj = 1;
    projarea = temp;
    //-- yz plane
    temp = (envelope.MaxY - envelope.MinY) * (envelope.MaxZ - envelope.MinZ); 
    if ( temp > projarea) {
      proj = 0;
    }
  }
  return proj;
}


bool triangulateandtag_xy(OGRGeometry* geometry, Triangulationxy &triangulation) {
  OGRPolygon *polygon = (OGRPolygon *)geometry;
  for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
    triangulation.insert_constraint(Point(polygon->getExteriorRing()->getX(currentPoint), 
                                          polygon->getExteriorRing()->getY(currentPoint),
                                          polygon->getExteriorRing()->getZ(currentPoint)),
                                    Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()), 
                                          polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                          polygon->getExteriorRing()->getZ((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
  } 
  for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
    for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint)
      triangulation.insert_constraint(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint), 
                                            polygon->getInteriorRing(currentRing)->getY(currentPoint),
                                            polygon->getInteriorRing(currentRing)->getZ(currentPoint)),
                                      Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()), 
                                            polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                            polygon->getInteriorRing(currentRing)->getZ((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
  }
  
  std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
  if (triangulation.number_of_faces() < 1) {
    return NULL;
  }
  
  // Tag
  void *interior = malloc(sizeof(void *));
  void *exterior = malloc(sizeof(void *));
  tag_xy(triangulation, interior, exterior);
  
  return true;
}
void tag_xy(Triangulationxy &triangulation, void *interiorHandle, void *exteriorHandle) {
	
  // Clean tags
  for (Triangulationxy::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
    currentFace->info() = NULL;
  
  // Initialise tagging
  std::stack<Triangulationxy::Face_handle> interiorStack, exteriorStack;
  exteriorStack.push(triangulation.infinite_face());
  std::stack<Triangulationxy::Face_handle> *currentStack = &exteriorStack;
  std::stack<Triangulationxy::Face_handle> *dualStack = &interiorStack;
  void *currentHandle = exteriorHandle;
  void *dualHandle = interiorHandle;
  
  // Until we finish
  while (!interiorStack.empty() || !exteriorStack.empty()) {
    
    // Give preference to whatever we're already doing
    while (!currentStack->empty()) {
      Triangulationxy::Face_handle currentFace = currentStack->top();
			currentStack->pop();
      if (currentFace->info() != NULL) continue;
			currentFace->info() = currentHandle;
      for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
        if (currentFace->neighbor(currentEdge)->info() == NULL)
          if (currentFace->is_constrained(currentEdge)) dualStack->push(currentFace->neighbor(currentEdge));
          else currentStack->push(currentFace->neighbor(currentEdge));
      }
    }
    
    // Flip
    if (currentHandle == exteriorHandle) {
      currentHandle = interiorHandle;
      dualHandle = exteriorHandle;
      currentStack = &interiorStack;
      dualStack = &exteriorStack;
    } else {
      currentHandle = exteriorHandle;
      dualHandle = interiorHandle;
      currentStack = &exteriorStack;
      dualStack = &interiorStack;
    }
	}
}

///////////////////////


bool triangulateandtag_xz(OGRGeometry* geometry, Triangulationxz &triangulation) {
  OGRPolygon *polygon = (OGRPolygon *)geometry;
  for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
    triangulation.insert_constraint(Point(polygon->getExteriorRing()->getX(currentPoint), 
                                          polygon->getExteriorRing()->getY(currentPoint),
                                          polygon->getExteriorRing()->getZ(currentPoint)),
                                    Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()), 
                                          polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                          polygon->getExteriorRing()->getZ((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
  } 
  for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
    for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint)
      triangulation.insert_constraint(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint), 
                                            polygon->getInteriorRing(currentRing)->getY(currentPoint),
                                            polygon->getInteriorRing(currentRing)->getZ(currentPoint)),
                                      Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()), 
                                            polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                            polygon->getInteriorRing(currentRing)->getZ((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
  }
  
  std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
  if (triangulation.number_of_faces() < 1) {
    return NULL;
  }
  
  // Tag
  void *interior = malloc(sizeof(void *));
  void *exterior = malloc(sizeof(void *));
  tag_xz(triangulation, interior, exterior);
  
  return true;
}
void tag_xz(Triangulationxz &triangulation, void *interiorHandle, void *exteriorHandle) {
	
  // Clean tags
  for (Triangulationxz::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
    currentFace->info() = NULL;
  
  // Initialise tagging
  std::stack<Triangulationxz::Face_handle> interiorStack, exteriorStack;
  exteriorStack.push(triangulation.infinite_face());
  std::stack<Triangulationxz::Face_handle> *currentStack = &exteriorStack;
  std::stack<Triangulationxz::Face_handle> *dualStack = &interiorStack;
  void *currentHandle = exteriorHandle;
  void *dualHandle = interiorHandle;
  
  // Until we finish
  while (!interiorStack.empty() || !exteriorStack.empty()) {
    
    // Give preference to whatever we're already doing
    while (!currentStack->empty()) {
      Triangulationxz::Face_handle currentFace = currentStack->top();
			currentStack->pop();
      if (currentFace->info() != NULL) continue;
			currentFace->info() = currentHandle;
      for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
        if (currentFace->neighbor(currentEdge)->info() == NULL)
          if (currentFace->is_constrained(currentEdge)) dualStack->push(currentFace->neighbor(currentEdge));
          else currentStack->push(currentFace->neighbor(currentEdge));
      }
    }
    
    // Flip
    if (currentHandle == exteriorHandle) {
      currentHandle = interiorHandle;
      dualHandle = exteriorHandle;
      currentStack = &interiorStack;
      dualStack = &exteriorStack;
    } else {
      currentHandle = exteriorHandle;
      dualHandle = interiorHandle;
      currentStack = &exteriorStack;
      dualStack = &interiorStack;
    }
	}
}


///////////////////////


bool triangulateandtag_yz(OGRGeometry* geometry, Triangulationyz &triangulation) {
  OGRPolygon *polygon = (OGRPolygon *)geometry;
  for (int currentPoint = 0; currentPoint < polygon->getExteriorRing()->getNumPoints(); ++currentPoint) {
    triangulation.insert_constraint(Point(polygon->getExteriorRing()->getX(currentPoint), 
                                          polygon->getExteriorRing()->getY(currentPoint),
                                          polygon->getExteriorRing()->getZ(currentPoint)),
                                    Point(polygon->getExteriorRing()->getX((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()), 
                                          polygon->getExteriorRing()->getY((currentPoint+1)%polygon->getExteriorRing()->getNumPoints()),
                                          polygon->getExteriorRing()->getZ((currentPoint+1)%polygon->getExteriorRing()->getNumPoints())));
  } 
  for (int currentRing = 0; currentRing < polygon->getNumInteriorRings(); ++currentRing) {
    for (int currentPoint = 0; currentPoint < polygon->getInteriorRing(currentRing)->getNumPoints(); ++currentPoint)
      triangulation.insert_constraint(Point(polygon->getInteriorRing(currentRing)->getX(currentPoint), 
                                            polygon->getInteriorRing(currentRing)->getY(currentPoint),
                                            polygon->getInteriorRing(currentRing)->getZ(currentPoint)),
                                      Point(polygon->getInteriorRing(currentRing)->getX((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()), 
                                            polygon->getInteriorRing(currentRing)->getY((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints()),
                                            polygon->getInteriorRing(currentRing)->getZ((currentPoint+1)%polygon->getInteriorRing(currentRing)->getNumPoints())));
  }
  
  std::cout << "Triangulation: " << triangulation.number_of_faces() << " faces, " << triangulation.number_of_vertices() << " vertices." << std::endl;
  if (triangulation.number_of_faces() < 1) {
    return NULL;
  }
  
  // Tag
  void *interior = malloc(sizeof(void *));
  void *exterior = malloc(sizeof(void *));
  tag_yz(triangulation, interior, exterior);
  
  return true;
}
void tag_yz(Triangulationyz &triangulation, void *interiorHandle, void *exteriorHandle) {
	
  // Clean tags
  for (Triangulationyz::Face_handle currentFace = triangulation.all_faces_begin(); currentFace != triangulation.all_faces_end(); ++currentFace)
    currentFace->info() = NULL;
  
  // Initialise tagging
  std::stack<Triangulationyz::Face_handle> interiorStack, exteriorStack;
  exteriorStack.push(triangulation.infinite_face());
  std::stack<Triangulationyz::Face_handle> *currentStack = &exteriorStack;
  std::stack<Triangulationyz::Face_handle> *dualStack = &interiorStack;
  void *currentHandle = exteriorHandle;
  void *dualHandle = interiorHandle;
  
  // Until we finish
  while (!interiorStack.empty() || !exteriorStack.empty()) {
    
    // Give preference to whatever we're already doing
    while (!currentStack->empty()) {
      Triangulationyz::Face_handle currentFace = currentStack->top();
			currentStack->pop();
      if (currentFace->info() != NULL) continue;
			currentFace->info() = currentHandle;
      for (int currentEdge = 0; currentEdge < 3; ++currentEdge) {
        if (currentFace->neighbor(currentEdge)->info() == NULL)
          if (currentFace->is_constrained(currentEdge)) dualStack->push(currentFace->neighbor(currentEdge));
          else currentStack->push(currentFace->neighbor(currentEdge));
      }
    }
    
    // Flip
    if (currentHandle == exteriorHandle) {
      currentHandle = interiorHandle;
      dualHandle = exteriorHandle;
      currentStack = &interiorStack;
      dualStack = &exteriorStack;
    } else {
      currentHandle = exteriorHandle;
      dualHandle = interiorHandle;
      currentStack = &exteriorStack;
      dualStack = &interiorStack;
    }
	}
}





