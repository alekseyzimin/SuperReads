#include <gtest/gtest.h>
#include <src/charb.hpp>

TEST(CharbBasic, Init) {
  charb empty;
  EXPECT_EQ((size_t)0, empty.capacity());
  EXPECT_EQ((size_t)0, empty.len());
  
  charb init_len(20);
  EXPECT_EQ((size_t)20, init_len.capacity());
  EXPECT_EQ((size_t)0, init_len.len());
  EXPECT_EQ('\0', init_len[0]);
  EXPECT_EQ('\0', *init_len);
  init_len.ensure(10);
  EXPECT_EQ((size_t)20, init_len.capacity());
  init_len.ensure(30);
  EXPECT_EQ((size_t)40, init_len.capacity());

  charb copy_const(init_len);
  EXPECT_EQ(init_len.capacity(), copy_const.capacity());
  EXPECT_EQ('\0', copy_const[0]);
  charb copy_op;
  copy_op = init_len;
  EXPECT_EQ(init_len.capacity(), copy_op.capacity());
  EXPECT_EQ('\0', copy_op[0]);
}

TEST(CharbBasic, Copy) {
  static const char *str = "Hello the world";
  const size_t str_len = strlen(str);
  const std::string string(str);

  charb from_str(str);
  EXPECT_EQ(str_len + 1, from_str.capacity());
  EXPECT_EQ(str_len, from_str.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_str[i]);

  charb from_str_len(str, str_len);
  EXPECT_EQ(str_len + 1, from_str_len.capacity());
  EXPECT_EQ(str_len, from_str_len.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_str_len[i]);

  charb from_string(string);
  EXPECT_EQ(str_len + 1, from_string.capacity());
  EXPECT_EQ(str_len, from_string.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_string[i]);
}

TEST(CharbBasic, Cast) {
  static const char *str = "This is a char buffer";
  const size_t str_len = strlen(str);
  charb b(str);
  
  char *s = b;
  EXPECT_EQ((void*)s, (void*)&b[0]);
  EXPECT_EQ(str_len, strlen(b));
  EXPECT_EQ(0, strcmp(str, b));
}

TEST(CharbBasic, fgets) {
  FILE *tf = tmpfile();
  if(!tf)
    throw std::runtime_error("Can't create tmp file");
  const char *l1 = "Hello there\n";
  const char *l2 = "A very long and meaningless sentence-line.";
  fputs(l1, tf); 
  fputs(l2, tf);
  rewind(tf);

  charb b(5);
  char *res = fgets(b, 3, tf); // 3 ignored!
  EXPECT_EQ(0, strcmp(l1, res));
  EXPECT_EQ(0, strcmp(l1, b));
  EXPECT_EQ(20, b.capacity());
  EXPECT_EQ(strlen(l1), b.len());

  res = fgets(b, tf);
  EXPECT_EQ(0, strcmp(l2, res));
  EXPECT_EQ(0, strcmp(l2, b));
  EXPECT_EQ(80, b.capacity());
  EXPECT_EQ(strlen(l2), b.len());

  res = fgets(b, tf);
  EXPECT_EQ((char *)0, res);
  EXPECT_EQ(0, strcmp(l2, b));
}
