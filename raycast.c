#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PLANE 0
#define SPHERE 1
#define CAMERA 2
#define LIGHT 3

#define MAX_COLOR_VALUE 255

int line = 1;

typedef struct {
  unsigned char r, g, b;
} Pixel;

typedef struct {
  double width;
  double height;
} Camera;

typedef struct {
  int kind; // 0 = Plane, 1 = Sphere, 2 = Camera, 3 = Light
  double diffuseColor[3];
  double specularColor[3];
  double position[3];
  union {
    struct {
      double normal[3];
    } plane;
    struct {
      double radius;
    } sphere;
  };
} Object;

typedef struct {
  double color[3];
  double position[3];
  double direction[3];
  double radialAtten[3];
  double angularAtten;
  double theta;
} Light;

Pixel* pixmap;
Camera** camera;
Object** objects;
Light** lights;

static inline double degreesToRads(double d) {
  return d * 0.0174533;
}

static inline double sqr(double v) {
  return v*v;
}

static inline void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

static inline double dot(const double* a, const double* b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double magnitude(const double* v) {
  return sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
}

static inline void subtract(double* v1, const double* v2) {
  v1[0] -= v2[0];
  v1[1] -= v2[1];
  v1[2] -= v2[2];
}

double clamp(double value, double min, double max) {
  if (value < min) return min;
  if (value > max) return max;
  return value;
}

void scale(double*v, double s) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

void reflect(const double* v, const double* n, double* r) {
  double dotResult = dot(v, n);
  dotResult *= 2;
  double nNew[3] = {
    n[0],
    n[1],
    n[2]
  };
  scale(nNew, dotResult);
  r[0] = v[0];
  r[1] = v[1];
  r[2] = v[2];
  subtract(r, nNew);
}

void negate(double* v) {
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}


// Wraps the getc() function and provides error checking and
// number maintenance
int fnextc(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("fnextc: '%c'\n", c);
#endif
  if (c == '\n') {
    line++;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}

// fexpectc() checks that the next character in d. If it is not it
// emits and error.
void fexpectc(FILE* json, int d) {
  int c = fnextc(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);
}

// skipWhitespace skips white space in the file.
void skipWhitespace(FILE* json) {
  int c = fnextc(json);
  while (isspace(c)) {
    c = fnextc(json);
  }
  ungetc(c, json);
}

// parseString gets the next string from the file handle and emits
// an error if a string can not be obtained.
char* nextString(FILE* json) {
  char buffer[129];
  int c = fnextc(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }
  c = fnextc(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported. See line %d.\n", line);
      exit(1);
    } else if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supperted. See line %d.\n", line);
      exit(1);
    } else if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters. See line %d.\n", line);
      exit(1);
    }
    buffer[i] = c;
    i++;
    c = fnextc(json);
  }
  buffer[i] = '\0';
  return strdup(buffer);
}

double nextNumber(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  // Error check this...
  return value;
}

double* nextVector(FILE* json) {
  double* v = malloc(3 * sizeof(double));
  fexpectc(json, '[');
  skipWhitespace(json);
  v[0] = nextNumber(json);
  skipWhitespace(json);
  fexpectc(json, ',');
  skipWhitespace(json);
  v[1] = nextNumber(json);
  skipWhitespace(json);
  fexpectc(json, ',');
  skipWhitespace(json);
  v[2] = nextNumber(json);
  skipWhitespace(json);
  fexpectc(json, ']');
  return v;
}

void parseObject(FILE* json, int currentObject, int objectType) {
  int c;

  if (objectType == SPHERE || objectType == PLANE) {
    objects[currentObject]->specularColor[0] = 0;
    objects[currentObject]->specularColor[1] = 0;
    objects[currentObject]->specularColor[2] = 0;
  }

  if (objectType == LIGHT) {
    lights[currentObject]->direction[0] = 0;
    lights[currentObject]->direction[1] = 0;
    lights[currentObject]->direction[2] = 0;
    lights[currentObject]->radialAtten[0] = INFINITY;
    lights[currentObject]->radialAtten[1] = INFINITY;
    lights[currentObject]->radialAtten[2] = INFINITY;
    lights[currentObject]->angularAtten = INFINITY;
    lights[currentObject]->theta = 0;
  }

  while (1) {
    c = fnextc(json);
    if (c == '}') {
      // Stop parsing this object

      if (objectType == LIGHT) {
        if (lights[currentObject]->radialAtten[0] != INFINITY || lights[currentObject]->radialAtten[1] != INFINITY || lights[currentObject]->radialAtten[2] != INFINITY) {
          if (lights[currentObject]->radialAtten[0] == INFINITY) {
            lights[currentObject]->radialAtten[0] = 0;
          }
          if (lights[currentObject]->radialAtten[1] == INFINITY) {
            lights[currentObject]->radialAtten[1] = 0;
          }
          if (lights[currentObject]->radialAtten[2] == INFINITY) {
            lights[currentObject]->radialAtten[2] = 1;
          }
        }
      }
      break;
    } else if (c == ',') {
      skipWhitespace(json);
      char* key = nextString(json);
      skipWhitespace(json);
      fexpectc(json, ':');
      skipWhitespace(json);

      if (strcmp(key, "width") == 0) {
        if (objectType == CAMERA) {
          double w = nextNumber(json);
          if (w > 0) {
              camera[0]->width = w;
          } else {
            fprintf(stderr, "Camera width must be greater than 0.\n");
            exit(1);
          }
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "height") == 0) {
        if (objectType == CAMERA) {
          double h = nextNumber(json);
          if (h > 0) {
              camera[0]->height = h;
          } else {
            fprintf(stderr, "Camera height must be greater than 0.\n");
            exit(1);
          }
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "radius") == 0) {
        if (objectType == SPHERE) {
          double radius = nextNumber(json);
          if (radius >= 0) {
            objects[currentObject]->sphere.radius = radius;
          } else {
            fprintf(stderr, "Error: Radius cannot be less than 0.\n");
            exit(1);
          }
        }  else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "color") == 0) {
        if (objectType == LIGHT) {
          double* v = nextVector(json);
          for (int i = 0; i < 3; i++) {
            lights[currentObject]->color[i] = v[i];
          }
          free(v);
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "diffuse_color") == 0) {
        if (objectType == PLANE || objectType == SPHERE) {
          double* v = nextVector(json);
          for (int i = 0; i < 3; i++) {
            objects[currentObject]->diffuseColor[i] = v[i];
          }
          free(v);
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      }  else if (strcmp(key, "specular_color") == 0) {
        if (objectType == PLANE || objectType == SPHERE) {
          double* v = nextVector(json);
          for (int i = 0; i < 3; i++) {
            objects[currentObject]->specularColor[i] = v[i];
          }
          free(v);
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "position") == 0) {
        if (objectType == PLANE || objectType == SPHERE) {
          double* v = nextVector(json);
          for (int i = 0; i < 3; i++) {
            objects[currentObject]->position[i] = v[i];
          }
          free(v);
        } else if (objectType == LIGHT) {
          double* v = nextVector(json);
          for (int i = 0; i < 3; i++) {
            lights[currentObject]->position[i] = v[i];
          }
          free(v);
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "normal") == 0) {
        if (objectType == PLANE) {
          double* v = nextVector(json);
          normalize(v);
          for (int i = 0; i < 3; i++) {
            objects[currentObject]->plane.normal[i] = v[i];
          }
          free(v);
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "direction") == 0) {
        if (objectType == LIGHT) {
          double* v = nextVector(json);
          normalize(v);
          for (int i = 0; i < 3; i++) {
            lights[currentObject]->direction[i] = v[i];
          }
          free(v);
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "radial-a2") == 0) {
        if (objectType == LIGHT) {
          double rAtten2 = nextNumber(json);
          lights[currentObject]->radialAtten[2] = rAtten2;
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "radial-a1") == 0) {
        if (objectType == LIGHT) {
          double rAtten1 = nextNumber(json);
          lights[currentObject]->radialAtten[1] = rAtten1;
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "radial-a0") == 0) {
        if (objectType == LIGHT) {
          double rAtten0 = nextNumber(json);
          lights[currentObject]->radialAtten[0] = rAtten0;
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "angular-a0") == 0) {
        if (objectType == LIGHT) {
          double aAtten = nextNumber(json);
          lights[currentObject]->angularAtten = aAtten;
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else if (strcmp(key, "theta") == 0) {
        if (objectType == LIGHT) {
          double theta = nextNumber(json);
          lights[currentObject]->theta = theta;
        } else {
          fprintf(stderr, "Error: Improper object field on line %d", line);
          exit(1);
        }
      } else {
        fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n", key, line);
        exit(1);
      }
      skipWhitespace(json);
    } else {
      fprintf(stderr, "Error: Unexpected value on line %d.\n", line);
      exit(1);
    }
  }
}

void parseJSON(char* fileName) {
  int c;
  FILE* json = fopen(fileName, "r");
  camera[0] = NULL;

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", fileName);
    exit(1);
  }

  skipWhitespace(json);

  // Find the beginning of the list
  fexpectc(json, '[');

  skipWhitespace(json);

  int currentObject = 0;
  int currentLight = 0;
  while (1) {
    c = fnextc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    } else if (c == '{') {
      skipWhitespace(json);

      // Parse the object
      char* key = nextString(json);
      if (strcmp(key, "type") != 0) {
        fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
        exit(1);
      }

      skipWhitespace(json);

      fexpectc(json, ':');

      skipWhitespace(json);

      char* value = nextString(json);

      skipWhitespace(json);
      if (strcmp(value, "camera") == 0) {
        if (camera[0] == NULL) {
          camera[0] = malloc(sizeof(Camera));
          parseObject(json, currentObject, CAMERA);
        } else {
          fprintf(stderr, "Error: There should only be one camera per scene.\n");
          exit(1);
        }
      } else if (strcmp(value, "sphere") == 0) {
        objects[currentObject] = malloc(sizeof(Object));
        objects[currentObject]->kind = SPHERE;
        parseObject(json, currentObject, SPHERE);
        currentObject++;
      } else if (strcmp(value, "plane") == 0) {
        objects[currentObject] = malloc(sizeof(Object));
        objects[currentObject]->kind = PLANE;
        parseObject(json, currentObject, PLANE);
        currentObject++;
      } else if (strcmp(value, "light") == 0) {
        lights[currentLight] = malloc(sizeof(Light));
        parseObject(json, currentLight, LIGHT);
        currentLight++;
      } else {
        fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
        exit(1);
      }

      skipWhitespace(json);
      c = fnextc(json);
      if (c == ',') {
        skipWhitespace(json);
      } else if (c == ']') {
        if (camera[0] == NULL) {
          fprintf(stderr, "Error: Scene must contain a camera.\n");
          exit(1);
        }
        objects[currentObject] = NULL;
        lights[currentLight] = NULL;
        fclose(json);
        return;
      } else {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        exit(1);
      }
    } else {
      fprintf(stderr, "Error: Expecting '{' on line %d.\n", line);
      exit(1);
    }
  }
}

double planeIntersection(const double* Ro, const double* Rd, const double* P, const double* N) {
  double Vd = dot(N, Rd);
  if (Vd == 0) return -1;
  double dist[3] = {
    P[0] - Ro[0],
    P[1] - Ro[1],
    P[2] - Ro[2]
  };
  double Vo = dot(dist, N);
  double t = Vo / Vd;
  if (t < 0) return -2;
  return t;
}

double sphereIntersection(const double* Ro, const double* Rd, const double* P, double r) {
  double A = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
  double B = 2 * (Rd[0] * (Ro[0] - P[0]) + Rd[1] * (Ro[1] - P[1]) + Rd[2] * (Ro[2] - P[2]));
  double C = sqr(Ro[0] - P[0]) + sqr(Ro[1] - P[1]) + sqr(Ro[2] - P[2]) - sqr(r);

  double det = sqr(B) - 4 * A * C;
  if (det < 0) return -1;
  det = sqrt(det);
  double t0 = (-B - det) / (2 * A);
  if (t0 > 0) return t0;

  double t1 = (-B + det) / (2 * A);
  if (t1 > 0) return t1;

  return -1;
}

double angularAttenuation(const double* Vo, const double* Vl, double a1, double angle) {
  double dotResult = dot(Vo, Vl);
  if (acos(dotResult) > angle / 2) {
    return 0;
  } else {
    return pow(dotResult, a1);
  }
}

double radialAttenuation(double a2, double a1, double a0, double d) {
  double quotient = a2 * sqr(d) + a1 * d + a0;
  if (quotient == 0) {
    return 0;
  }
  if (d == INFINITY) {
    return 1;
  } else {
    return 1.0 / quotient;
  }
}

double diffuseReflection(double Kd, double Il, const double* N, const double* L) {
  double dotResult = dot(N, L);
  if (dotResult > 0) {
    return Kd * Il * dotResult;
  } else {
    return 0;
  }
}

double specularReflection(double Ks, double Il, const double* V, const double* R, const double* N, const double* L, double ns) {
  double dotResult = dot(V, R);
  if (dotResult > 0 && dot(N, L) > 0) {
    return Ks * Il * pow(dotResult, ns);
  } else {
    return 0;
  }
}

void createScene(int width, int height) {
  double cx = 0;
  double cy = 0;
  double h = camera[0]->height;
  double w = camera[0]->width;

  int M = height;
  int N = width;

  double pixheight = h / M;
  double pixwidth = w / N;

  for (int y = 0; y < M; y++) {
    for (int x = 0; x < N; x++) {
      double Ro[3] = {0, 0, 0};
      double Rd[3] = {
        cx - (w/2) + pixwidth * (x + 0.5),
        cy - (h/2) + pixheight * (y + 0.5),
        1
      };
      normalize(Rd);

      double closestT = INFINITY;
      Object* closestObject = NULL;
      for (int i = 0; objects[i] != NULL; i++) {
        double t = 0;

        switch(objects[i]->kind) {
          case PLANE:
            t = planeIntersection(Ro, Rd,
              objects[i]->position,
              objects[i]->plane.normal);
            break;
          case SPHERE:
            t = sphereIntersection(Ro, Rd,
              objects[i]->position,
              objects[i]->sphere.radius);
            break;
          default:
            fprintf(stderr, "Error: Object does not have an appropriate kind.");
            exit(1);
        }

        if (t > 0 && t < closestT) {
          closestT = t;
          closestObject = objects[i];
        }
      }

      if (closestT < INFINITY) {
        double color[3];
        color[0] = 0;
        color[1] = 0;
        color[2] = 0;

        for (int i = 0; lights[i] != NULL; i++) {
          double RoNew[3] = {
            closestT * Rd[0] + Ro[0],
            closestT * Rd[1] + Ro[1],
            closestT * Rd[2] + Ro[2]
          };
          double RdNew[3] = {
            lights[i]->position[0] - RoNew[0],
            lights[i]->position[1] - RoNew[1],
            lights[i]->position[2] - RoNew[2]
          };

          normalize(RdNew);

          int shadow = 0;
          for (int j = 0; objects[j] != NULL; j++) {
            double t = 0;
            if (objects[j] == closestObject) continue;
            switch(objects[j]->kind) {
              case PLANE:
                t = planeIntersection(RoNew, RdNew,
                  objects[j]->position,
                  objects[j]->plane.normal);
                break;
              case SPHERE:
                t = sphereIntersection(RoNew, RdNew,
                  objects[j]->position,
                  objects[j]->sphere.radius);
                break;
              default:
                fprintf(stderr, "Error: Object does not have an appropriate kind.");
                exit(1);
            }
            if (t > 0 && t < magnitude(RdNew)) {
              shadow = 1;
              break;
            }
          }

          if (shadow == 0) {
            double N[3];
            if (closestObject->kind == PLANE) {
              N[0] = closestObject->plane.normal[0];
              N[1] = closestObject->plane.normal[1];
              N[2] = closestObject->plane.normal[2];
            } else if (closestObject->kind == SPHERE) {
              N[0] = RoNew[0] - closestObject->position[0];
              N[1] = RoNew[1] - closestObject->position[1];
              N[2] = RoNew[2] - closestObject->position[2];
            }

            normalize(N);
            double L[3] = {
              RdNew[0],
              RdNew[1],
              RdNew[2]
            };
            normalize(L);
            double LNeg[3] = {
              -L[0],
              -L[1],
              -L[2]
            };
            double R[3];
            reflect(L, N, R);
            double V[3] = {
              Rd[0],
              Rd[1],
              Rd[2]
            };

            double pos[3] = {
              lights[i]->position[0],
              lights[i]->position[1],
              lights[i]->position[2]
            };
            subtract(pos, RoNew);
            double d = magnitude(pos);

            double col;
            for (int c = 0; c < 3; c++) {
              col = 1;
              if (lights[i]->angularAtten != INFINITY && lights[i]->theta != 0) {
                col *= angularAttenuation(LNeg, lights[i]->direction, lights[i]->angularAtten, degreesToRads(lights[i]->theta));
              }
              if (lights[i]->radialAtten[0] != INFINITY) {
                col *= radialAttenuation(lights[i]->radialAtten[2], lights[i]->radialAtten[1], lights[i]->radialAtten[0], d);
              }
              col *= (diffuseReflection(closestObject->diffuseColor[c], lights[i]->color[c], N, L) + (specularReflection(closestObject->specularColor[c], lights[i]->color[c], V, R, N, L, 20)));
              color[c] += col;
            }
          }
        }
        if (closestObject != NULL) {
          pixmap[(M - 1) * N - (y * N) + x].r = (unsigned char)(clamp(color[0], 0, 1) * MAX_COLOR_VALUE);
          pixmap[(M - 1) * N - (y * N) + x].g = (unsigned char)(clamp(color[1], 0, 1) * MAX_COLOR_VALUE);
          pixmap[(M - 1) * N - (y * N) + x].b = (unsigned char)(clamp(color[2], 0, 1) * MAX_COLOR_VALUE);
        }  else {
            pixmap[(M - 1) * N - (y * N) + x].r = 0;
            pixmap[(M - 1) * N - (y * N) + x].g = 0;
            pixmap[(M - 1) * N - (y * N) + x].b = 0;
          }
    } else {
        pixmap[(M - 1) * N - (y * N) + x].r = 0;
        pixmap[(M - 1) * N - (y * N) + x].g = 0;
        pixmap[(M - 1) * N - (y * N) + x].b = 0;
      }
    }
  }
}

void writeP6(char* outputPath, int width, int height) {
  FILE* fh = fopen(outputPath, "wb");
  if (fh == NULL) {
    fprintf(stderr, "Error: Output file not found.\n");
  }
  fprintf(fh, "P6\n# Converted with Robert Rasmussen's ppmrw\n%d %d\n%d\n", width, height, MAX_COLOR_VALUE);
  fwrite(pixmap, sizeof(Pixel), width*height, fh);
  fclose(fh);
}

int main(int argc, char* argv[]) {
  if (argc != 5) {
    fprintf(stderr, "Usage: raycast width height input.json output.ppm");
    exit(1);
  }

  int width = atoi(argv[1]);
  if (width <= 0) {
    fprintf(stderr, "Error: Width must be greater than 0.");
    exit(1);
  }
  int height = atoi(argv[2]);
  if (height <= 0) {
    fprintf(stderr, "Error: Width must be greater than 0.");
    exit(1);
  }

  pixmap = malloc(sizeof(Pixel) * width * height);
  camera = malloc(sizeof(Camera));
  objects = malloc(sizeof(Object*) * 129);
  lights = malloc(sizeof(Light*) * 129);

  parseJSON(argv[3]);
  createScene(width, height);

  writeP6(argv[4], width, height);

  free(camera[0]);
  free(camera);

  for (int i = 0; objects[i] != NULL; i++) {
    free(objects[i]);
  }
  free(objects);

  for (int i = 0; lights[i] != NULL; i++) {
    free(lights[i]);
  }
  free(lights);

#ifdef DEBUG
  displayObjects();
#endif


  return 0;
}
