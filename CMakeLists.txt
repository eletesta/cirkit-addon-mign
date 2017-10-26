add_cirkit_library(
  NAME cirkit_mign
  AUTO_DIRS src
  USE
    cirkit_classical
    cirkit_cli
    cirkit_formal_z3
  INCLUDE
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
  DEFINE
    PUBLIC ADDON_MIGN
  OPTIONS
    PRIVATE -Wno-sign-compare
  COMMANDS
    classical/cli/commands/mign_commands.hpp
)

if( enable_cirkit-addon-formal )
  set(FORMAL_LIBRARY cirkit_formal_z3)
else()
  set(FORMAL_LIBRARY)
endif()