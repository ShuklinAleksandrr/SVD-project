/**
 * @brief Defines an elliptical region and checks if a point is inside.
 * 
 * This program uses SLEPc to define an elliptical region and checks whether
 * the point (0.1, 0.3) is inside or outside that region.
 */


static char help[] = "Simple example.\n\n";

#include <slepcrg.h>

int main(int argc, char **argv) {
  PetscErrorCode ierr;
  RG rg;
  PetscInt inside;
  PetscReal re, im;
  PetscScalar ar, ai;

  ierr = SlepcInitialize(&argc, &argv, (char*)0, help);
  if (ierr) return ierr;

  ierr = RGCreate(PETSC_COMM_WORLD, &rg);
  CHKERRQ(ierr);
  ierr = RGSetType(rg, RGELLIPSE);
  CHKERRQ(ierr);
  ierr = RGEllipseSetParameters(rg, 1.1, 2, 0.1);
  CHKERRQ(ierr);
  ierr = RGSetFromOptions(rg);
  CHKERRQ(ierr);
  
  re = 0.1;
  im = 0.3;

#ifdef PETSC_USE_COMPLEX
  ar = re + im * PETSC_i;
#else
  ar = re;
  ai = im;
#endif

  ierr = RGCheckInside(rg, 1, &ar, &ai, &inside);
  CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Point (%g, %g) is %s the region\n", (double)re, (double)im, inside ? "inside" : "outside");
  ierr = RGDestroy(&rg); 
  CHKERRQ(ierr);

  ierr = SlepcFinalize();
  
  return ierr;
}