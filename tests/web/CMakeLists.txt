#SPDX-License-Identifier: LGPL-3.0-or-later
## Test web

# add the tests
ADD_TEST(web.smoke web.sh)
SET_TESTS_PROPERTIES(web
                     PROPERTIES LABELS smoke)
