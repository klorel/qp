#include "NlReader.h"
#include "ConvexRelaxation.h"

int main(int argc, char**argv) {
  std::string const fileName(argv[1]);
  Problem problem;
  NlReader nlReader(fileName, "");
  nlReader.createObjects(problem);
  std::cout << problem << std::endl;
//	std::map<std::size_t, double> result;
//	csdpInterface.bestReformulation(result);
  ConvexRelaxation convexRelaxation(problem);
  convexRelaxation.identifyBlocks();
  return 0;
}
