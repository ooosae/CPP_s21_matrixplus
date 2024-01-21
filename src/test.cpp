#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

#define EPS 1e-3

TEST(Test, test_base_constructor_1) {
  S21Matrix my_matrix;
  ASSERT_EQ(my_matrix.GetCols(), 0);
  ASSERT_EQ(my_matrix.GetRows(), 0);
}

TEST(Test, test_param_constructor_1) {
  S21Matrix my_matrix(2, 2);
  ASSERT_EQ(my_matrix.GetCols(), 2);
  ASSERT_EQ(my_matrix.GetRows(), 2);
  for (int i = 0; i < my_matrix.GetRows(); i++)
    for (int j = 0; j < my_matrix.GetCols(); j++)
      ASSERT_EQ((my_matrix(i, j)), 0);
}

TEST(Test, test_param_constructor_2) {
  S21Matrix my_matrix(10, 10);
  ASSERT_EQ(my_matrix.GetCols(), 10);
  ASSERT_EQ(my_matrix.GetRows(), 10);
}

TEST(Test, test_param_constructor_3) {
  S21Matrix my_matrix(2, 5);
  ASSERT_EQ(my_matrix.GetCols(), 5);
  ASSERT_EQ(my_matrix.GetRows(), 2);
}

TEST(Test, test_param_constructor_4) {
  ASSERT_ANY_THROW(S21Matrix my_matrix(-3, -2));
}

TEST(Test, test_param_constructor_5) {
  ASSERT_ANY_THROW(S21Matrix my_matrix(0, 0));
}

TEST(Test, test_param_constructor_6) {
  ASSERT_ANY_THROW(S21Matrix my_matrix(0, 1));
}

TEST(Test, test_param_constructor_7) {
  ASSERT_ANY_THROW(S21Matrix my_matrix(5, -5));
}

TEST(Test, test_param_constructor_8) {
  ASSERT_NO_THROW(S21Matrix my_matrix(5, 5));
}

//

TEST(Test, test_copy_constructor_1) {
  S21Matrix my_matrix(5, 5);
  S21Matrix new_matrix(my_matrix);

  ASSERT_EQ(new_matrix.GetCols(), my_matrix.GetCols());
  ASSERT_EQ(new_matrix.GetRows(), my_matrix.GetCols());
}

TEST(Test, test_copy_constructor_2) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 1;
  my_matrix(0, 2) = 1;
  my_matrix(1, 0) = 2;
  my_matrix(1, 1) = 2;
  my_matrix(1, 2) = 2;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix = my_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(new_matrix(i, j), my_matrix(i, j));
    }
  }
}

TEST(Test, test_copy_constructor_3) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.14159261;
  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;

  S21Matrix new_matrix = my_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(new_matrix(i, j), my_matrix(i, j));
    }
  }
}
TEST(Test, test_copy_constructor_4) {
  S21Matrix my_matrix(3, 5);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.1415926;
  my_matrix(0, 3) = 3.1415926;
  my_matrix(0, 4) = 3.1415926;

  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(1, 3) = 3.1415926;
  my_matrix(1, 4) = 3.1415926;

  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;
  my_matrix(2, 3) = 3.1415926;
  my_matrix(2, 4) = 3.1415926;

  S21Matrix new_matrix = my_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 5; j++) {
      ASSERT_EQ(new_matrix(i, j), my_matrix(i, j));
    }
  }
}

TEST(Test, test_copy_constructor_5) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.1415926;
  my_matrix(0, 3) = 3.1415926;
  my_matrix(0, 4) = 3.1415926;

  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(1, 3) = 3.1415926;
  my_matrix(1, 4) = 3.1415926;

  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;
  my_matrix(2, 3) = 3.1415926;
  my_matrix(2, 4) = 3.1415926;

  my_matrix(3, 0) = 3.1415926;
  my_matrix(3, 1) = 3.1415926;
  my_matrix(3, 2) = 3.1415926;
  my_matrix(3, 3) = 3.1415926;
  my_matrix(3, 4) = 3.1415926;

  my_matrix(4, 0) = 3.1415926;
  my_matrix(4, 1) = 3.1415926;
  my_matrix(4, 2) = 3.1415926;
  my_matrix(4, 3) = 3.1415926;
  my_matrix(4, 4) = 3.1415926;

  S21Matrix new_matrix = my_matrix;

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      ASSERT_EQ(new_matrix(i, j), my_matrix(i, j));
    }
  }
}

TEST(Test, test_replace_constructor_1) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 1;
  my_matrix(0, 2) = 1;
  my_matrix(1, 0) = 1;
  my_matrix(1, 1) = 1;
  my_matrix(1, 2) = 1;
  my_matrix(2, 0) = 1;
  my_matrix(2, 1) = 1;
  my_matrix(2, 2) = 1;

  S21Matrix new_matrix(std::move(my_matrix));

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(new_matrix(i, j), 1);
    }
  }
}

TEST(Test, test_replace_constructor_2) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.1415926;
  my_matrix(0, 3) = 3.1415926;
  my_matrix(0, 4) = 3.1415926;

  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(1, 3) = 3.1415926;
  my_matrix(1, 4) = 3.1415926;

  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;
  my_matrix(2, 3) = 3.1415926;
  my_matrix(2, 4) = 3.1415926;

  my_matrix(3, 0) = 3.1415926;
  my_matrix(3, 1) = 3.1415926;
  my_matrix(3, 2) = 3.1415926;
  my_matrix(3, 3) = 3.1415926;
  my_matrix(3, 4) = 3.1415926;

  my_matrix(4, 0) = 3.1415926;
  my_matrix(4, 1) = 3.1415926;
  my_matrix(4, 2) = 3.1415926;
  my_matrix(4, 3) = 3.1415926;
  my_matrix(4, 4) = 3.1415926;

  S21Matrix new_matrix(std::move(my_matrix));

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      ASSERT_EQ(new_matrix(i, j), 3.1415926);
    }
  }
}

TEST(Test, test_getters_setters) {
  S21Matrix my_matrix(3, 3);

  ASSERT_EQ(my_matrix.GetRows(), 3);
  ASSERT_EQ(my_matrix.GetCols(), 3);

  ASSERT_ANY_THROW(my_matrix.SetRows(0));
  ASSERT_ANY_THROW(my_matrix.SetRows(-1));
  ASSERT_ANY_THROW(my_matrix.SetCols(0));
  ASSERT_ANY_THROW(my_matrix.SetCols(-5));

  my_matrix.SetCols(5);

  ASSERT_EQ(my_matrix.GetRows(), 3);
  ASSERT_EQ(my_matrix.GetCols(), 5);
}
/*
TEST(Test, test_getters_setters_1) {
    S21Matrix my_matrix(3, 3);

    ASSERT_EQ(my_matrix.GetRows(), 3);
    ASSERT_EQ(my_matrix.GetCols(), 3);

    my_matrix.SetRows(6);
    my_matrix.SetCols(5);

    ASSERT_EQ(my_matrix.GetRows(), 6);
    ASSERT_EQ(my_matrix.GetCols(), 5);

    my_matrix.SetRows(2);
    my_matrix.SetCols(2);

    ASSERT_EQ(my_matrix.GetRows(), 2);
    ASSERT_EQ(my_matrix.GetCols(), 2);

    my_matrix.SetRows(10);
    my_matrix.SetCols(10);

    ASSERT_EQ(my_matrix.GetRows(), 10);
    ASSERT_EQ(my_matrix.GetCols(), 10);
}
*/
TEST(Test, test_eq_matrix_1) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.14159261;
  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = 3.1415926;
  new_matrix(0, 1) = 3.1415926;
  new_matrix(0, 2) = 3.14159261;
  new_matrix(1, 0) = 3.1415926;
  new_matrix(1, 1) = 3.1415926;
  new_matrix(1, 2) = 3.1415926;
  new_matrix(2, 0) = 3.1415926;
  new_matrix(2, 1) = 3.1415926;
  new_matrix(2, 2) = 3.1415926;

  bool answer = (new_matrix == my_matrix);
  ASSERT_TRUE(answer);
}

TEST(Test, test_eq_matrix_2) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.14159261;
  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;

  S21Matrix new_matrix(5, 5);
  new_matrix(0, 0) = 3.1415926;
  new_matrix(0, 1) = 3.1415926;
  new_matrix(0, 2) = 3.14159261;
  new_matrix(1, 0) = 3.1415926;
  new_matrix(1, 1) = 3.1415926;
  new_matrix(1, 2) = 3.1415926;
  new_matrix(2, 0) = 3.1415926;
  new_matrix(2, 1) = 3.1415926;
  new_matrix(2, 2) = 3.1415926;
  new_matrix(3, 0) = 3.1415926;
  new_matrix(3, 1) = 3.1415926;
  new_matrix(3, 2) = 3.1415926;
  new_matrix(3, 3) = 3.1415926;
  new_matrix(3, 4) = 3.1415926;
  new_matrix(4, 0) = 3.1415926;
  new_matrix(4, 1) = 3.1415926;
  new_matrix(4, 2) = 3.1415926;
  new_matrix(4, 3) = 3.1415926;
  new_matrix(4, 4) = 3.1415926;

  bool answer = (new_matrix == my_matrix);
  ASSERT_FALSE(answer);
}

TEST(Test, test_eq_matrix_3) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.14159261;
  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;

  S21Matrix new_matrix(3, 2);
  new_matrix(0, 0) = 3.1415926;
  new_matrix(0, 1) = 3.1415926;
  new_matrix(1, 0) = 3.1415926;
  new_matrix(1, 1) = 3.1415926;
  new_matrix(2, 0) = 3.1415926;
  new_matrix(2, 1) = 3.1415926;

  bool answer = (new_matrix == my_matrix);
  ASSERT_FALSE(answer);
}

TEST(Test, test_eq_matrix_4) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.14159261;
  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = 3.1415926;
  new_matrix(0, 1) = 3.1415926;
  new_matrix(0, 2) = 3.14159261;
  new_matrix(1, 0) = 3.1415926;
  new_matrix(1, 1) = 3.1415926;
  new_matrix(1, 2) = 3.1415926;
  new_matrix(2, 0) = 3.1415926;
  new_matrix(2, 1) = 3.1415926;
  new_matrix(2, 2) = 3.1415926;

  bool answer = (new_matrix == my_matrix);
  ASSERT_FALSE(answer);
}

TEST(Test, test_eq_matrix_5) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.14159261;
  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = 3.1415926;
  new_matrix(0, 1) = 5;
  new_matrix(0, 2) = 3.14159261;
  new_matrix(1, 0) = 3.1415926;
  new_matrix(1, 1) = 3.1415926;
  new_matrix(1, 2) = 3.1415926;
  new_matrix(2, 0) = 3.1415926;
  new_matrix(2, 1) = 0;
  new_matrix(2, 2) = 3.1415926;

  bool answer = (new_matrix == my_matrix);
  ASSERT_FALSE(answer);
}

TEST(Test, test_eq_matrix_6) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 3;
  my_matrix(1, 2) = 3;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = 3;
  new_matrix(0, 1) = 3;
  new_matrix(0, 2) = 3;
  new_matrix(1, 0) = 3;
  new_matrix(1, 1) = 3;
  new_matrix(1, 2) = 3;
  new_matrix(2, 0) = 3;
  new_matrix(2, 1) = 3;
  new_matrix(2, 2) = 3;

  bool answer = (new_matrix.EqMatrix(my_matrix));
  ASSERT_TRUE(answer);
}

TEST(Test, test_sum_matrix_1) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 3;
  my_matrix(1, 2) = 3;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = 3;
  new_matrix(0, 1) = 3;
  new_matrix(0, 2) = 3;
  new_matrix(1, 0) = 3;
  new_matrix(1, 1) = 3;
  new_matrix(1, 2) = 3;
  new_matrix(2, 0) = 3;
  new_matrix(2, 1) = 3;
  new_matrix(2, 2) = 3;

  my_matrix += new_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), 6);
    }
  }
}

TEST(Test, test_sum_matrix_2) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 3;
  my_matrix(1, 2) = 3;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  my_matrix += new_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), -2);
    }
  }
}

TEST(Test, test_sum_matrix_3) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 5;
  my_matrix(1, 2) = 6;
  my_matrix(2, 0) = 7;
  my_matrix(2, 1) = 8;
  my_matrix(2, 2) = 9;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  my_matrix += new_matrix;

  ASSERT_EQ(my_matrix(0, 0), -4);
  ASSERT_EQ(my_matrix(0, 1), -3);
  ASSERT_EQ(my_matrix(0, 2), -2);
  ASSERT_EQ(my_matrix(1, 0), -1);
  ASSERT_EQ(my_matrix(1, 1), 0);
  ASSERT_EQ(my_matrix(1, 2), 1);
  ASSERT_EQ(my_matrix(2, 0), 2);
  ASSERT_EQ(my_matrix(2, 1), 3);
  ASSERT_EQ(my_matrix(2, 2), 4);
}

TEST(Test, test_sum_matrix_4) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 5;
  my_matrix(1, 2) = 6;
  my_matrix(2, 0) = 7;
  my_matrix(2, 1) = 8;
  my_matrix(2, 2) = 9;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  S21Matrix res = my_matrix + new_matrix;

  ASSERT_EQ(res(0, 0), -4);
  ASSERT_EQ(res(0, 1), -3);
  ASSERT_EQ(res(0, 2), -2);
  ASSERT_EQ(res(1, 0), -1);
  ASSERT_EQ(res(1, 1), 0);
  ASSERT_EQ(res(1, 2), 1);
  ASSERT_EQ(res(2, 0), 2);
  ASSERT_EQ(res(2, 1), 3);
  ASSERT_EQ(res(2, 2), 4);
}

TEST(Test, test_sum_matrix_5) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 3;
  my_matrix(1, 2) = 3;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  my_matrix += new_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), -2);
    }
  }
}

TEST(Test, test_sum_matrix_6) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.1415926;
  my_matrix(0, 3) = 3.1415926;
  my_matrix(0, 4) = 3.1415926;

  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(1, 3) = 3.1415926;
  my_matrix(1, 4) = 3.1415926;

  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;
  my_matrix(2, 3) = 3.1415926;
  my_matrix(2, 4) = 3.1415926;

  my_matrix(3, 0) = 3.1415926;
  my_matrix(3, 1) = 3.1415926;
  my_matrix(3, 2) = 3.1415926;
  my_matrix(3, 3) = 3.1415926;
  my_matrix(3, 4) = 3.1415926;

  my_matrix(4, 0) = 3.1415926;
  my_matrix(4, 1) = 3.1415926;
  my_matrix(4, 2) = 3.1415926;
  my_matrix(4, 3) = 3.1415926;
  my_matrix(4, 4) = 3.1415926;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  ASSERT_ANY_THROW(my_matrix += new_matrix);
}

TEST(Test, test_sum_matrix_7) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.1415926;
  my_matrix(0, 3) = 3.1415926;
  my_matrix(0, 4) = 3.1415926;

  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(1, 3) = 3.1415926;
  my_matrix(1, 4) = 3.1415926;

  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;
  my_matrix(2, 3) = 3.1415926;
  my_matrix(2, 4) = 3.1415926;

  my_matrix(3, 0) = 3.1415926;
  my_matrix(3, 1) = 3.1415926;
  my_matrix(3, 2) = 3.1415926;
  my_matrix(3, 3) = 3.1415926;
  my_matrix(3, 4) = 3.1415926;

  my_matrix(4, 0) = 3.1415926;
  my_matrix(4, 1) = 3.1415926;
  my_matrix(4, 2) = 3.1415926;
  my_matrix(4, 3) = 3.1415926;
  my_matrix(4, 4) = 3.1415926;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  ASSERT_ANY_THROW((my_matrix.SumMatrix(new_matrix)));
}

TEST(Test, test_sub_matrix_1) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 3;
  my_matrix(1, 2) = 3;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = 3;
  new_matrix(0, 1) = 3;
  new_matrix(0, 2) = 3;
  new_matrix(1, 0) = 3;
  new_matrix(1, 1) = 3;
  new_matrix(1, 2) = 3;
  new_matrix(2, 0) = 3;
  new_matrix(2, 1) = 3;
  new_matrix(2, 2) = 3;

  my_matrix -= new_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), 0);
    }
  }
}

TEST(Test, test_sub_matrix_2) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 3;
  my_matrix(1, 2) = 3;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  my_matrix -= new_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), 8);
    }
  }
}

TEST(Test, test_sub_matrix_3) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 5;
  my_matrix(1, 2) = 6;
  my_matrix(2, 0) = 7;
  my_matrix(2, 1) = 8;
  my_matrix(2, 2) = 9;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  my_matrix -= new_matrix;

  ASSERT_EQ(my_matrix(0, 0), 6);
  ASSERT_EQ(my_matrix(0, 1), 7);
  ASSERT_EQ(my_matrix(0, 2), 8);
  ASSERT_EQ(my_matrix(1, 0), 9);
  ASSERT_EQ(my_matrix(1, 1), 10);
  ASSERT_EQ(my_matrix(1, 2), 11);
  ASSERT_EQ(my_matrix(2, 0), 12);
  ASSERT_EQ(my_matrix(2, 1), 13);
  ASSERT_EQ(my_matrix(2, 2), 14);
}

TEST(Test, test_sub_matrix_4) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 5;
  my_matrix(1, 2) = 6;
  my_matrix(2, 0) = 7;
  my_matrix(2, 1) = 8;
  my_matrix(2, 2) = 9;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  S21Matrix res = my_matrix - new_matrix;

  ASSERT_EQ(res(0, 0), 6);
  ASSERT_EQ(res(0, 1), 7);
  ASSERT_EQ(res(0, 2), 8);
  ASSERT_EQ(res(1, 0), 9);
  ASSERT_EQ(res(1, 1), 10);
  ASSERT_EQ(res(1, 2), 11);
  ASSERT_EQ(res(2, 0), 12);
  ASSERT_EQ(res(2, 1), 13);
  ASSERT_EQ(res(2, 2), 14);
}

TEST(Test, test_sub_matrix_5) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 3;
  my_matrix(0, 2) = 3;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 3;
  my_matrix(1, 2) = 3;
  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 3;
  my_matrix(2, 2) = 3;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  my_matrix -= new_matrix;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), 8);
    }
  }
}

TEST(Test, test_sub_matrix_6) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.1415926;
  my_matrix(0, 3) = 3.1415926;
  my_matrix(0, 4) = 3.1415926;

  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(1, 3) = 3.1415926;
  my_matrix(1, 4) = 3.1415926;

  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;
  my_matrix(2, 3) = 3.1415926;
  my_matrix(2, 4) = 3.1415926;

  my_matrix(3, 0) = 3.1415926;
  my_matrix(3, 1) = 3.1415926;
  my_matrix(3, 2) = 3.1415926;
  my_matrix(3, 3) = 3.1415926;
  my_matrix(3, 4) = 3.1415926;

  my_matrix(4, 0) = 3.1415926;
  my_matrix(4, 1) = 3.1415926;
  my_matrix(4, 2) = 3.1415926;
  my_matrix(4, 3) = 3.1415926;
  my_matrix(4, 4) = 3.1415926;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  ASSERT_ANY_THROW(my_matrix -= new_matrix);
}

TEST(Test, test_sub_matrix_7) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 3.1415926;
  my_matrix(0, 1) = 3.1415926;
  my_matrix(0, 2) = 3.1415926;
  my_matrix(0, 3) = 3.1415926;
  my_matrix(0, 4) = 3.1415926;

  my_matrix(1, 0) = 3.1415926;
  my_matrix(1, 1) = 3.1415926;
  my_matrix(1, 2) = 3.1415926;
  my_matrix(1, 3) = 3.1415926;
  my_matrix(1, 4) = 3.1415926;

  my_matrix(2, 0) = 3.1415926;
  my_matrix(2, 1) = 3.1415926;
  my_matrix(2, 2) = 3.1415926;
  my_matrix(2, 3) = 3.1415926;
  my_matrix(2, 4) = 3.1415926;

  my_matrix(3, 0) = 3.1415926;
  my_matrix(3, 1) = 3.1415926;
  my_matrix(3, 2) = 3.1415926;
  my_matrix(3, 3) = 3.1415926;
  my_matrix(3, 4) = 3.1415926;

  my_matrix(4, 0) = 3.1415926;
  my_matrix(4, 1) = 3.1415926;
  my_matrix(4, 2) = 3.1415926;
  my_matrix(4, 3) = 3.1415926;
  my_matrix(4, 4) = 3.1415926;

  S21Matrix new_matrix(3, 3);
  new_matrix(0, 0) = -5;
  new_matrix(0, 1) = -5;
  new_matrix(0, 2) = -5;
  new_matrix(1, 0) = -5;
  new_matrix(1, 1) = -5;
  new_matrix(1, 2) = -5;
  new_matrix(2, 0) = -5;
  new_matrix(2, 1) = -5;
  new_matrix(2, 2) = -5;

  ASSERT_ANY_THROW((my_matrix.SubMatrix(new_matrix)));
}

TEST(Test, test_mul_number_1) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 1;
  my_matrix(0, 2) = 1;
  my_matrix(1, 0) = 1;
  my_matrix(1, 1) = 1;
  my_matrix(1, 2) = 1;
  my_matrix(2, 0) = 1;
  my_matrix(2, 1) = 1;
  my_matrix(2, 2) = 1;

  S21Matrix result = my_matrix * 2;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(result(i, j), 2);
    }
  }
}

TEST(Test, test_mul_number_2) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 1;
  my_matrix(0, 2) = 1;
  my_matrix(1, 0) = 1;
  my_matrix(1, 1) = 1;
  my_matrix(1, 2) = 1;
  my_matrix(2, 0) = 1;
  my_matrix(2, 1) = 1;
  my_matrix(2, 2) = 1;

  S21Matrix result = my_matrix * -5;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(result(i, j), -5);
    }
  }
}

TEST(Test, test_mul_number_3) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 0.1;
  my_matrix(0, 1) = 0.1;
  my_matrix(0, 2) = 0.1;
  my_matrix(1, 0) = 0.1;
  my_matrix(1, 1) = 0.1;
  my_matrix(1, 2) = 0.1;
  my_matrix(2, 0) = 0.1;
  my_matrix(2, 1) = 0.1;
  my_matrix(2, 2) = 0.1;

  S21Matrix result = my_matrix * -5;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(result(i, j), -0.5);
    }
  }
}

TEST(Test, test_mul_number_4) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 0.1;
  my_matrix(0, 1) = 0.1;
  my_matrix(0, 2) = 0.1;
  my_matrix(1, 0) = 0.1;
  my_matrix(1, 1) = 0.1;
  my_matrix(1, 2) = 0.1;
  my_matrix(2, 0) = 0.1;
  my_matrix(2, 1) = 0.1;
  my_matrix(2, 2) = 0.1;

  my_matrix *= -5;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), -0.5);
    }
  }
}

TEST(Test, test_mul_number_5) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 0.5;
  my_matrix(0, 1) = 0.5;
  my_matrix(0, 2) = 0.5;
  my_matrix(1, 0) = 0.5;
  my_matrix(1, 1) = 0.5;
  my_matrix(1, 2) = 0.5;
  my_matrix(2, 0) = 0.5;
  my_matrix(2, 1) = 0.5;
  my_matrix(2, 2) = 0.5;

  my_matrix *= -5;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), -2.5);
    }
  }
}

TEST(Test, test_mul_number_6) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 0.5;
  my_matrix(0, 1) = 0.5;
  my_matrix(0, 2) = 0.5;
  my_matrix(1, 0) = 0.5;
  my_matrix(1, 1) = 0.5;
  my_matrix(1, 2) = 0.5;
  my_matrix(2, 0) = 0.5;
  my_matrix(2, 1) = 0.5;
  my_matrix(2, 2) = 0.5;

  my_matrix *= 5;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), 2.5);
    }
  }
}

TEST(Test, test_mul_number_7) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 0.5;
  my_matrix(0, 1) = 0.5;
  my_matrix(0, 2) = 0.5;
  my_matrix(1, 0) = 0.5;
  my_matrix(1, 1) = 0.5;
  my_matrix(1, 2) = 0.5;
  my_matrix(2, 0) = 0.5;
  my_matrix(2, 1) = 0.5;
  my_matrix(2, 2) = 0.5;

  S21Matrix result(3, 3);

  result = my_matrix * -5;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(result(i, j), -2.5);
    }
  }
}

TEST(Test, test_mul_number_8) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 0.5;
  my_matrix(0, 1) = 0.5;
  my_matrix(0, 2) = 0.5;
  my_matrix(1, 0) = 0.5;
  my_matrix(1, 1) = 0.5;
  my_matrix(1, 2) = 0.5;
  my_matrix(2, 0) = 0.5;
  my_matrix(2, 1) = 0.5;
  my_matrix(2, 2) = 0.5;

  my_matrix.MulNumber(5);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(my_matrix(i, j), 2.5);
    }
  }
}

TEST(Test, test_mul_matrix_1) {
  S21Matrix my_matrix(3, 2);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 4;

  my_matrix(1, 0) = 2;
  my_matrix(1, 1) = 5;

  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 6;

  S21Matrix new_matrix(2, 3);
  new_matrix(0, 0) = 1;
  new_matrix(0, 1) = -1;
  new_matrix(0, 2) = 1;

  new_matrix(1, 0) = 2;
  new_matrix(1, 1) = 3;
  new_matrix(1, 2) = 4;

  my_matrix *= new_matrix;

  ASSERT_EQ(my_matrix(0, 0), 9);
  ASSERT_EQ(my_matrix(0, 1), 11);
  ASSERT_EQ(my_matrix(0, 2), 17);
  ASSERT_EQ(my_matrix(1, 0), 12);
  ASSERT_EQ(my_matrix(1, 1), 13);
  ASSERT_EQ(my_matrix(1, 2), 22);
  ASSERT_EQ(my_matrix(2, 0), 15);
  ASSERT_EQ(my_matrix(2, 1), 15);
  ASSERT_EQ(my_matrix(2, 2), 27);
}

TEST(Test, test_mul_matrix_2) {
  S21Matrix my_matrix(3, 1);
  S21Matrix new_matrix(2, 2);

  ASSERT_ANY_THROW((my_matrix *= new_matrix));
}

TEST(Test, test_mul_matrix_3) {
  S21Matrix my_matrix(3, 1);
  my_matrix(0, 0) = 1.0;
  my_matrix(1, 0) = 2.0;
  my_matrix(2, 0) = 3.0;

  S21Matrix new_matrix(1, 3);
  new_matrix(0, 0) = 9.0;
  new_matrix(0, 1) = 11.0;
  new_matrix(0, 2) = 17.0;

  my_matrix *= new_matrix;

  ASSERT_EQ(my_matrix(0, 0), 9.0);
  ASSERT_EQ(my_matrix(0, 1), 11.0);
  ASSERT_EQ(my_matrix(0, 2), 17.0);
  ASSERT_EQ(my_matrix(1, 0), 18.0);
  ASSERT_EQ(my_matrix(1, 1), 22.0);
  ASSERT_EQ(my_matrix(1, 2), 34.0);
  ASSERT_EQ(my_matrix(2, 0), 27.0);
  ASSERT_EQ(my_matrix(2, 1), 33.0);
  ASSERT_EQ(my_matrix(2, 2), 51.0);
}

TEST(Test, test_mul_matrix_4) {
  S21Matrix my_matrix(3, 1);
  my_matrix(0, 0) = 1.0;
  my_matrix(1, 0) = 2.0;
  my_matrix(2, 0) = 3.0;

  S21Matrix new_matrix(1, 3);
  new_matrix(0, 0) = 9.0;
  new_matrix(0, 1) = 11.0;
  new_matrix(0, 2) = 17.0;

  my_matrix = my_matrix * new_matrix;

  ASSERT_EQ(my_matrix(0, 0), 9.0);
  ASSERT_EQ(my_matrix(0, 1), 11.0);
  ASSERT_EQ(my_matrix(0, 2), 17.0);
  ASSERT_EQ(my_matrix(1, 0), 18.0);
  ASSERT_EQ(my_matrix(1, 1), 22.0);
  ASSERT_EQ(my_matrix(1, 2), 34.0);
  ASSERT_EQ(my_matrix(2, 0), 27.0);
  ASSERT_EQ(my_matrix(2, 1), 33.0);
  ASSERT_EQ(my_matrix(2, 2), 51.0);
}

TEST(Test, test_mul_matrix_5) {
  S21Matrix my_matrix(3, 2);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 4;

  my_matrix(1, 0) = 2;
  my_matrix(1, 1) = 5;

  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 6;

  S21Matrix new_matrix(2, 3);
  new_matrix(0, 0) = 1;
  new_matrix(0, 1) = -1;
  new_matrix(0, 2) = 1;

  new_matrix(1, 0) = 2;
  new_matrix(1, 1) = 3;
  new_matrix(1, 2) = 4;

  my_matrix.MulMatrix(new_matrix);

  ASSERT_EQ(my_matrix(0, 0), 9);
  ASSERT_EQ(my_matrix(0, 1), 11);
  ASSERT_EQ(my_matrix(0, 2), 17);
  ASSERT_EQ(my_matrix(1, 0), 12);
  ASSERT_EQ(my_matrix(1, 1), 13);
  ASSERT_EQ(my_matrix(1, 2), 22);
  ASSERT_EQ(my_matrix(2, 0), 15);
  ASSERT_EQ(my_matrix(2, 1), 15);
  ASSERT_EQ(my_matrix(2, 2), 27);
}

TEST(Test, test_mul_matrix_6) {
  ASSERT_ANY_THROW(S21Matrix my_matrix(3, 0));
  S21Matrix my_matrix(1, 1);
  S21Matrix new_matrix(2, 2);

  ASSERT_ANY_THROW((my_matrix *= new_matrix));
}

TEST(Test, test_trans_matrix_1) {
  S21Matrix my_matrix(3, 2);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 4;

  my_matrix(1, 0) = 2;
  my_matrix(1, 1) = 5;

  my_matrix(2, 0) = 3;
  my_matrix(2, 1) = 6;

  S21Matrix result = my_matrix.Transpose();

  ASSERT_EQ(result(0, 0), 1);
  ASSERT_EQ(result(0, 1), 2);
  ASSERT_EQ(result(0, 2), 3);
  ASSERT_EQ(result(1, 0), 4);
  ASSERT_EQ(result(1, 1), 5);
  ASSERT_EQ(result(1, 2), 6);
}

TEST(Test, test_trans_matrix_2) {
  S21Matrix my_matrix(5, 2);
  my_matrix(0, 0) = 3;
  my_matrix(0, 1) = 1;

  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 2;

  my_matrix(2, 0) = 5;
  my_matrix(2, 1) = 6;

  my_matrix(3, 0) = 3;
  my_matrix(3, 1) = 4;

  my_matrix(4, 0) = 8;
  my_matrix(4, 1) = 9;

  S21Matrix result = my_matrix.Transpose();

  ASSERT_EQ(result(0, 0), 3);
  ASSERT_EQ(result(0, 1), 4);
  ASSERT_EQ(result(0, 2), 5);
  ASSERT_EQ(result(0, 3), 3);
  ASSERT_EQ(result(0, 4), 8);
  ASSERT_EQ(result(1, 0), 1);
  ASSERT_EQ(result(1, 1), 2);
  ASSERT_EQ(result(1, 2), 6);
  ASSERT_EQ(result(1, 3), 4);
  ASSERT_EQ(result(1, 4), 9);
}

TEST(Test, test_trans_matrix_3) {
  S21Matrix my_matrix(3, 2);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 4;
  my_matrix(2, 0) = 5;
  my_matrix(2, 1) = 6;

  S21Matrix result = my_matrix.Transpose();

  ASSERT_EQ(result(0, 0), 1);
  ASSERT_EQ(result(0, 1), 3);
  ASSERT_EQ(result(0, 2), 5);
  ASSERT_EQ(result(1, 0), 2);
  ASSERT_EQ(result(1, 1), 4);
  ASSERT_EQ(result(1, 2), 6);
}

TEST(Test, test_trans_matrix_4) {
  S21Matrix my_matrix(2, 2);
  my_matrix(0, 0) = 4;
  my_matrix(0, 1) = 2;

  my_matrix(1, 0) = 9;
  my_matrix(1, 1) = 0;

  S21Matrix result = my_matrix.Transpose();

  ASSERT_EQ(result(0, 0), 4);
  ASSERT_EQ(result(0, 1), 9);

  ASSERT_EQ(result(1, 0), 2);
  ASSERT_EQ(result(1, 1), 0);
}

TEST(Test, test_trans_matrix_5) {
  S21Matrix my_matrix(3, 6);
  my_matrix(0, 0) = 5;
  my_matrix(0, 1) = 4;
  my_matrix(0, 2) = 5;
  my_matrix(0, 3) = 4;
  my_matrix(0, 4) = 3;
  my_matrix(0, 5) = 2;

  my_matrix(1, 0) = 1;
  my_matrix(1, 1) = 2;
  my_matrix(1, 2) = 3;
  my_matrix(1, 3) = 4;
  my_matrix(1, 4) = 5;
  my_matrix(1, 5) = 6;

  my_matrix(2, 0) = 7;
  my_matrix(2, 1) = 8;
  my_matrix(2, 2) = 9;
  my_matrix(2, 3) = 0;
  my_matrix(2, 4) = 1;
  my_matrix(2, 5) = 2;

  S21Matrix result = my_matrix.Transpose();

  ASSERT_EQ(result(0, 0), 5);
  ASSERT_EQ(result(0, 1), 1);
  ASSERT_EQ(result(0, 2), 7);

  ASSERT_EQ(result(1, 0), 4);
  ASSERT_EQ(result(1, 1), 2);
  ASSERT_EQ(result(1, 2), 8);

  ASSERT_EQ(result(2, 0), 5);
  ASSERT_EQ(result(2, 1), 3);
  ASSERT_EQ(result(2, 2), 9);

  ASSERT_EQ(result(3, 0), 4);
  ASSERT_EQ(result(3, 1), 4);
  ASSERT_EQ(result(3, 2), 0);

  ASSERT_EQ(result(4, 0), 3);
  ASSERT_EQ(result(4, 1), 5);
  ASSERT_EQ(result(4, 2), 1);

  ASSERT_EQ(result(5, 0), 2);
  ASSERT_EQ(result(5, 1), 6);
  ASSERT_EQ(result(5, 2), 2);
}

TEST(Test, test_determinant_matrix_1) {
  S21Matrix my_matrix(4, 4);
  my_matrix(0, 0) = 9;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 2;
  my_matrix(0, 3) = 4;

  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 4;
  my_matrix(1, 2) = 4;
  my_matrix(1, 3) = 4;

  my_matrix(2, 0) = 4;
  my_matrix(2, 1) = 4;
  my_matrix(2, 2) = 9;
  my_matrix(2, 3) = 9;

  my_matrix(3, 0) = 1;
  my_matrix(3, 1) = 1;
  my_matrix(3, 2) = 5;
  my_matrix(3, 3) = 1;

  double det = my_matrix.Determinant();

  ASSERT_EQ(det, -578);
}

TEST(Test, test_determinant_matrix_2) {
  S21Matrix my_matrix(2, 2);
  my_matrix(0, 0) = 4;
  my_matrix(0, 1) = 1;

  my_matrix(1, 0) = 2;
  my_matrix(1, 1) = 5;

  double det = my_matrix.Determinant();

  ASSERT_EQ(det, 18);
}

TEST(Test, test_determinant_matrix_3) {
  S21Matrix my_matrix(2, 2);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 1;

  my_matrix(1, 0) = 1;
  my_matrix(1, 1) = 1;

  double det = my_matrix.Determinant();

  ASSERT_EQ(det, 0);
}

TEST(Test, test_determinant_matrix_4) {
  S21Matrix my_matrix(2, 2);
  my_matrix(0, 0) = 5;
  my_matrix(0, 1) = 1;

  my_matrix(1, 0) = 1;
  my_matrix(1, 1) = 5;

  double det = my_matrix.Determinant();

  ASSERT_EQ(det, 24);
}

TEST(Test, test_determinant_matrix_5) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 1;
  my_matrix(0, 3) = 3;
  my_matrix(0, 4) = 1;

  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 1;
  my_matrix(1, 2) = 5;
  my_matrix(1, 3) = 1;
  my_matrix(1, 4) = 6;

  my_matrix(2, 0) = 1;
  my_matrix(2, 1) = 7;
  my_matrix(2, 2) = 1;
  my_matrix(2, 3) = 8;
  my_matrix(2, 4) = 1;

  my_matrix(3, 0) = 9;
  my_matrix(3, 1) = 1;
  my_matrix(3, 2) = 0;
  my_matrix(3, 3) = 1;
  my_matrix(3, 4) = 2;

  my_matrix(4, 0) = 2;
  my_matrix(4, 1) = 2;
  my_matrix(4, 2) = 3;
  my_matrix(4, 3) = 2;
  my_matrix(4, 4) = 4;

  int det = my_matrix.Determinant();

  ASSERT_DOUBLE_EQ(det, 110);
}

TEST(Test, test_determinant_matrix_6) {
  S21Matrix my_matrix(4, 5);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 1;
  my_matrix(0, 3) = 3;
  my_matrix(0, 4) = 1;

  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 1;
  my_matrix(1, 2) = 5;
  my_matrix(1, 3) = 1;
  my_matrix(1, 4) = 6;

  my_matrix(2, 0) = 1;
  my_matrix(2, 1) = 7;
  my_matrix(2, 2) = 1;
  my_matrix(2, 3) = 8;
  my_matrix(2, 4) = 1;

  my_matrix(3, 0) = 9;
  my_matrix(3, 1) = 1;
  my_matrix(3, 2) = 0;
  my_matrix(3, 3) = 1;
  my_matrix(3, 4) = 2;

  ASSERT_ANY_THROW(my_matrix.Determinant());
}

TEST(Test, test_determinant_matrix_7) {
  S21Matrix my_matrix(1, 1);
  my_matrix(0, 0) = 1;

  double det = my_matrix.Determinant();

  ASSERT_EQ(det, 1);
}

TEST(Test, test_determinant_matrix_8) {
  S21Matrix my_matrix(6, 6);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 1;
  my_matrix(0, 2) = 1;
  my_matrix(0, 3) = 1;
  my_matrix(0, 4) = 1;
  my_matrix(0, 5) = 3;

  my_matrix(1, 0) = 2;
  my_matrix(1, 1) = 2;
  my_matrix(1, 2) = 2;
  my_matrix(1, 3) = 2;
  my_matrix(1, 4) = 2;
  my_matrix(1, 5) = 3;

  my_matrix(2, 0) = 4;
  my_matrix(2, 1) = 4;
  my_matrix(2, 2) = 4;
  my_matrix(2, 3) = 4;
  my_matrix(2, 4) = 4;
  my_matrix(2, 5) = 5;

  my_matrix(3, 0) = 5;
  my_matrix(3, 1) = 5;
  my_matrix(3, 2) = 5;
  my_matrix(3, 3) = 5;
  my_matrix(3, 4) = 5;
  my_matrix(3, 5) = 6;

  my_matrix(4, 0) = 6;
  my_matrix(4, 1) = 6;
  my_matrix(4, 2) = 6;
  my_matrix(4, 3) = 6;
  my_matrix(4, 4) = 6;
  my_matrix(4, 5) = 7;

  my_matrix(5, 0) = 7;
  my_matrix(5, 1) = 7;
  my_matrix(5, 2) = 7;
  my_matrix(5, 3) = 7;
  my_matrix(5, 4) = 7;
  my_matrix(5, 5) = 8;

  double det = my_matrix.Determinant();

  ASSERT_EQ(det, 0);
}

TEST(Test, test_determinant_matrix_9) {
  S21Matrix my_matrix(6, 6);
  my_matrix(0, 0) = -11;
  my_matrix(0, 1) = 11;
  my_matrix(0, 2) = 12;
  my_matrix(0, 3) = 10;
  my_matrix(0, 4) = 12;
  my_matrix(0, 5) = 1;

  my_matrix(1, 0) = 31;
  my_matrix(1, 1) = 55;
  my_matrix(1, 2) = 44;
  my_matrix(1, 3) = 22;
  my_matrix(1, 4) = 15;
  my_matrix(1, 5) = 0;

  my_matrix(2, 0) = 4;
  my_matrix(2, 1) = 4;
  my_matrix(2, 2) = 4;
  my_matrix(2, 3) = 4;
  my_matrix(2, 4) = 4;
  my_matrix(2, 5) = 5;

  my_matrix(3, 0) = -33;
  my_matrix(3, 1) = 5;
  my_matrix(3, 2) = 5;
  my_matrix(3, 3) = 18;
  my_matrix(3, 4) = 5;
  my_matrix(3, 5) = 6;

  my_matrix(4, 0) = 6;
  my_matrix(4, 1) = 6;
  my_matrix(4, 2) = 42;
  my_matrix(4, 3) = 6;
  my_matrix(4, 4) = 6;
  my_matrix(4, 5) = 7;

  my_matrix(5, 0) = 7;
  my_matrix(5, 1) = 7;
  my_matrix(5, 2) = 7;
  my_matrix(5, 3) = 7;
  my_matrix(5, 4) = 7;
  my_matrix(5, 5) = 8;

  int det = my_matrix.Determinant();

  ASSERT_EQ(det, -1568808);
}

TEST(Test, test_calc_compliments_1) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 3;

  my_matrix(1, 0) = 0;
  my_matrix(1, 1) = 4;
  my_matrix(1, 2) = 2;

  my_matrix(2, 0) = 5;
  my_matrix(2, 1) = 2;
  my_matrix(2, 2) = 1;

  S21Matrix result = my_matrix.CalcComplements();

  ASSERT_EQ(result(0, 0), 0);
  ASSERT_EQ(result(0, 1), 10);
  ASSERT_EQ(result(0, 2), -20);
  ASSERT_EQ(result(1, 0), 4);
  ASSERT_EQ(result(1, 1), -14);
  ASSERT_EQ(result(1, 2), 8);
  ASSERT_EQ(result(2, 0), -8);
  ASSERT_EQ(result(2, 1), -2);
  ASSERT_EQ(result(2, 2), 4);
}

TEST(Test, test_calc_compliments_2) {
  S21Matrix my_matrix(4, 4);
  my_matrix(0, 0) = 9;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 2;
  my_matrix(0, 3) = 4;

  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 4;
  my_matrix(1, 2) = 4;
  my_matrix(1, 3) = 4;

  my_matrix(2, 0) = 4;
  my_matrix(2, 1) = 4;
  my_matrix(2, 2) = 9;
  my_matrix(2, 3) = 9;

  my_matrix(3, 0) = 1;
  my_matrix(3, 1) = 1;
  my_matrix(3, 2) = 5;
  my_matrix(3, 3) = 1;

  S21Matrix result;
  result = my_matrix.CalcComplements();

  ASSERT_NEAR(result(0, 0), -80, EPS);
  ASSERT_NEAR(result(0, 1), 44, EPS);
  ASSERT_NEAR(result(0, 2), 5, EPS);
  ASSERT_NEAR(result(0, 3), 11, EPS);

  ASSERT_NEAR(result(1, 0), 18, EPS);
  ASSERT_NEAR(result(1, 1), -270, EPS);
  ASSERT_NEAR(result(1, 2), 35, EPS);
  ASSERT_NEAR(result(1, 3), 77, EPS);

  ASSERT_NEAR(result(2, 0), 32, EPS);
  ASSERT_NEAR(result(2, 1), 98, EPS);
  ASSERT_NEAR(result(2, 2), -2, EPS);
  ASSERT_NEAR(result(2, 3), -120, EPS);

  ASSERT_NEAR(result(3, 0), -40, EPS);
  ASSERT_NEAR(result(3, 1), 22, EPS);
  ASSERT_NEAR(result(3, 2), -142, EPS);
  ASSERT_NEAR(result(3, 3), 150, EPS);
}

TEST(Test, test_calc_compliments_3) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = 9;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 2;
  my_matrix(0, 3) = 4;
  my_matrix(0, 4) = 1;

  my_matrix(1, 0) = 3;
  my_matrix(1, 1) = 4;
  my_matrix(1, 2) = 4;
  my_matrix(1, 3) = 4;
  my_matrix(1, 4) = 1;

  my_matrix(2, 0) = 4;
  my_matrix(2, 1) = 4;
  my_matrix(2, 2) = 9;
  my_matrix(2, 3) = 9;
  my_matrix(2, 4) = 5;

  my_matrix(3, 0) = 1;
  my_matrix(3, 1) = 1;
  my_matrix(3, 2) = 5;
  my_matrix(3, 3) = 1;
  my_matrix(3, 4) = 3;

  my_matrix(4, 0) = 1;
  my_matrix(4, 1) = 1;
  my_matrix(4, 2) = 5;
  my_matrix(4, 3) = 1;
  my_matrix(4, 4) = 2;

  S21Matrix result = my_matrix.CalcComplements();

  ASSERT_NEAR(result(0, 0), 80, EPS);
  ASSERT_NEAR(result(0, 1), -44, EPS);
  ASSERT_NEAR(result(0, 2), -5, EPS);
  ASSERT_NEAR(result(0, 3), -11, EPS);
  ASSERT_NEAR(result(0, 4), 0, EPS);

  ASSERT_NEAR(result(1, 0), -18, EPS);
  ASSERT_NEAR(result(1, 1), 270, EPS);
  ASSERT_NEAR(result(1, 2), -35, EPS);
  ASSERT_NEAR(result(1, 3), -77, EPS);
  ASSERT_NEAR(result(1, 4), -0, EPS);

  ASSERT_NEAR(result(2, 0), -32, EPS);
  ASSERT_NEAR(result(2, 1), -98, EPS);
  ASSERT_NEAR(result(2, 2), 2, EPS);
  ASSERT_NEAR(result(2, 3), 120, EPS);
  ASSERT_NEAR(result(2, 4), 0, EPS);

  ASSERT_NEAR(result(3, 0), 18, EPS);
  ASSERT_NEAR(result(3, 1), 308, EPS);
  ASSERT_NEAR(result(3, 2), -254, EPS);
  ASSERT_NEAR(result(3, 3), -212, EPS);
  ASSERT_NEAR(result(3, 4), 578, EPS);

  ASSERT_NEAR(result(4, 0), 22, EPS);
  ASSERT_NEAR(result(4, 1), -330, EPS);
  ASSERT_NEAR(result(4, 2), 396, EPS);
  ASSERT_NEAR(result(4, 3), 62, EPS);
  ASSERT_NEAR(result(4, 4), -578, EPS);
}

TEST(Test, test_calc_compliments_4) {
  S21Matrix my_matrix(1, 2);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 2;

  ASSERT_ANY_THROW(my_matrix.CalcComplements());
}

TEST(Test, test_inverse_matrix_1) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 4;
  my_matrix(0, 1) = -2;
  my_matrix(0, 2) = 1;

  my_matrix(1, 0) = 1;
  my_matrix(1, 1) = 6;
  my_matrix(1, 2) = -2;

  my_matrix(2, 0) = 1;
  my_matrix(2, 1) = 0;
  my_matrix(2, 2) = 0;

  S21Matrix for_eq(3, 3);

  for_eq(0, 0) = 1;
  for_eq(0, 1) = 0;
  for_eq(0, 2) = 0;

  for_eq(1, 0) = 0;
  for_eq(1, 1) = 1;
  for_eq(1, 2) = 0;

  for_eq(2, 0) = 0;
  for_eq(2, 1) = 0;
  for_eq(2, 2) = 1;

  S21Matrix new_matrix(3, 3);
  new_matrix = my_matrix;

  S21Matrix result = my_matrix.InverseMatrix();
  new_matrix *= result;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(new_matrix(i, j), for_eq(i, j));
    }
  }
}

TEST(Test, test_inverse_matrix_2) {
  S21Matrix my_matrix(3, 3);
  my_matrix(0, 0) = 8;
  my_matrix(0, 1) = 2;
  my_matrix(0, 2) = 4;

  my_matrix(1, 0) = 4;
  my_matrix(1, 1) = 6;
  my_matrix(1, 2) = 6;

  my_matrix(2, 0) = 4;
  my_matrix(2, 1) = 8;
  my_matrix(2, 2) = 8;

  S21Matrix for_eq(3, 3);

  for_eq(0, 0) = 1;
  for_eq(0, 1) = 0;
  for_eq(0, 2) = 0;

  for_eq(1, 0) = 0;
  for_eq(1, 1) = 1;
  for_eq(1, 2) = 0;

  for_eq(2, 0) = 0;
  for_eq(2, 1) = 0;
  for_eq(2, 2) = 1;

  S21Matrix new_matrix(3, 3);
  new_matrix = my_matrix;

  S21Matrix result = my_matrix.InverseMatrix();
  new_matrix *= result;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ASSERT_EQ(new_matrix(i, j), for_eq(i, j));
    }
  }
}

TEST(Test, test_inverse_matrix_3) {
  S21Matrix my_matrix(5, 5);
  my_matrix(0, 0) = -1.5;
  my_matrix(0, 1) = 2.2;
  my_matrix(0, 2) = 4.6;
  my_matrix(0, 3) = 2.2;
  my_matrix(0, 4) = -4.6;

  my_matrix(1, 0) = -4;
  my_matrix(1, 1) = -66;
  my_matrix(1, 2) = 3;
  my_matrix(1, 3) = 27.6;
  my_matrix(1, 4) = 3.2;

  my_matrix(2, 0) = 4.3;
  my_matrix(2, 1) = 8.9;
  my_matrix(2, 2) = -62.8;
  my_matrix(2, 3) = -100;
  my_matrix(2, 4) = 100;

  my_matrix(3, 0) = 14;
  my_matrix(3, 1) = 18;
  my_matrix(3, 2) = -6;
  my_matrix(3, 3) = -10;
  my_matrix(3, 4) = 1;

  my_matrix(4, 0) = 42.5;
  my_matrix(4, 1) = 81.9;
  my_matrix(4, 2) = -6.8;
  my_matrix(4, 3) = 1.5;
  my_matrix(4, 4) = 77;

  S21Matrix for_eq(5, 5);

  for_eq(0, 0) = 1;
  for_eq(0, 1) = 0;
  for_eq(0, 2) = 0;
  for_eq(0, 3) = 0;
  for_eq(0, 4) = 0;

  for_eq(1, 0) = 0;
  for_eq(1, 1) = 1;
  for_eq(1, 2) = 0;
  for_eq(1, 3) = 0;
  for_eq(1, 4) = 0;

  for_eq(2, 0) = 0;
  for_eq(2, 1) = 0;
  for_eq(2, 2) = 1;
  for_eq(2, 3) = 0;
  for_eq(2, 4) = 0;

  for_eq(3, 0) = 0;
  for_eq(3, 1) = 0;
  for_eq(3, 2) = 0;
  for_eq(3, 3) = 1;
  for_eq(3, 4) = 0;

  for_eq(4, 0) = 0;
  for_eq(4, 1) = 0;
  for_eq(4, 2) = 0;
  for_eq(4, 3) = 0;
  for_eq(4, 4) = 1;

  S21Matrix new_matrix(5, 5);
  new_matrix = my_matrix;

  S21Matrix result = my_matrix.InverseMatrix();
  new_matrix *= result;

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      ASSERT_NEAR(new_matrix(i, j), for_eq(i, j), 1e-7);
    }
  }
}

GTEST_API_ int main() {
  testing::InitGoogleTest();
  return RUN_ALL_TESTS();
}