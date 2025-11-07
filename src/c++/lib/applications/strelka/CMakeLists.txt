#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

include(${THIS_CXX_LIBRARY_CMAKE})

# Link OpenMP to the strelka library if available
# This is needed because strelka_run.cpp uses OpenMP
if (OpenMP_CXX_FOUND)
    target_link_libraries(${LIBRARY_TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(${LIBRARY_TARGET_NAME} PUBLIC USE_OPENMP=1)
    message(STATUS "Linking strelka library (${LIBRARY_TARGET_NAME}) with OpenMP")
    message(STATUS "Defining USE_OPENMP=1 for ${LIBRARY_TARGET_NAME}")
else()
    target_compile_definitions(${LIBRARY_TARGET_NAME} PUBLIC USE_OPENMP=0)
    message(STATUS "OpenMP not available, USE_OPENMP=0 for ${LIBRARY_TARGET_NAME}")
endif()
