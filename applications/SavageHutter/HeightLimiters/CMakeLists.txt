include_directories(.)

add_library(HeightLimiters STATIC
  PositiveLayerLimiter.cpp
  SqueezeLimiterWithLayer.cpp
  AverageValuesNonNegativeLimiter.cpp
)
