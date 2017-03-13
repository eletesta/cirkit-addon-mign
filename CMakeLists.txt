add_cirkit_library(
  NAME cirkit_mign
  AUTO_DIRS src
  USE
    cirkit_classical
  INCLUDE
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
  DEFINE
    PUBLIC ADDON_MIGN
  OPTIONS
    PRIVATE -Wno-sign-compare
  COMMANDS
    classical/cli/commands/mign_commands.hpp
)
