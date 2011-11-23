#include <gtest/gtest.h>
#include <src/charb.hpp>

TEST(CharbBasic, Init) {
  charb empty;
  EXPECT_EQ((size_t)0, empty.capacity());
  EXPECT_EQ((size_t)0, empty.len());
  
  charb init_len(20);
  EXPECT_EQ((size_t)21, init_len.capacity());
  EXPECT_EQ((size_t)0, init_len.len());
  EXPECT_EQ('\0', init_len[0]);
  EXPECT_EQ('\0', *init_len);
  init_len.reserve(10);
  EXPECT_EQ((size_t)21, init_len.capacity());
  init_len.reserve(30);
  EXPECT_EQ((size_t)42, init_len.capacity());

  charb copy_const(init_len);
  EXPECT_EQ(init_len.len(), copy_const.len());
  EXPECT_EQ('\0', copy_const[0]);
  charb copy_op;
  copy_op = init_len;
  EXPECT_EQ(init_len.len(), copy_op.len());
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
  EXPECT_STREQ(str, from_str);
  EXPECT_EQ('H', from_str.front());
  EXPECT_EQ('d', from_str.back());

  charb from_str_len(str, str_len);
  EXPECT_EQ(str_len + 1, from_str_len.capacity());
  EXPECT_EQ(str_len, from_str_len.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_str_len[i]);
  EXPECT_STREQ(str, from_str_len);
  EXPECT_EQ('H', from_str_len.front());
  EXPECT_EQ('d', from_str_len.back());

  charb from_string(string);
  EXPECT_EQ(str_len + 1, from_string.capacity());
  EXPECT_EQ(str_len, from_string.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_string[i]);
  EXPECT_STREQ(str, from_string);
  EXPECT_EQ('H', from_string.front());
  EXPECT_EQ('d', from_string.back());

  charb from_copy;
  from_copy = from_string;
  EXPECT_EQ(str_len + 1, from_copy.capacity());
  EXPECT_EQ(str_len, from_copy.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_copy[i]);
  EXPECT_STREQ(str, from_copy);
  EXPECT_EQ('H', from_copy.front());
  EXPECT_EQ('d', from_copy.back());

  charb from_copy_str;
  from_copy_str = str;
  EXPECT_EQ(str_len + 1, from_copy_str.capacity());
  EXPECT_EQ(str_len, from_copy_str.len());
  for(size_t i = 0; i < str_len; ++i)
    EXPECT_EQ(str[i], from_copy_str[i]);
  EXPECT_STREQ(str, from_copy_str);
  EXPECT_EQ('H', from_copy_str.front());
  EXPECT_EQ('d', from_copy_str.back());
  

  charb x("Hello");
  x[3] = '_';
  EXPECT_STREQ("Hel_o", x);
}

TEST(CharbBasic, Cast) {
  static const char *str = "This is a char buffer";
  const size_t str_len = strlen(str);
  charb b(str);
  
  char *s = b;
  EXPECT_EQ((void*)s, (void*)&b[0]);
  EXPECT_EQ(str_len, strlen(b));
  EXPECT_STREQ(str, b);
}

class CharbStd : public ::testing::Test {
public:
  virtual void SetUp() {
    s1 = "Hello you";
    s2 = "How are you";
  }
  const char *s1, *s2;
};

TEST_F(CharbStd, Strcat) {
  charb b;
  char *res;

  res = strcat(b, s1);
  EXPECT_STREQ(s1, b);
  EXPECT_EQ(strlen(s1), b.len());
  
  std::string sres(s1);
  sres += s2;
  res = strcat(b, s2);
  EXPECT_STREQ(sres.c_str(), b);
  EXPECT_EQ(sres.size(), b.len());
}

TEST_F(CharbStd, Strcpy) {
  charb b;
  char *res;

  res = strcpy(b, s1);
  EXPECT_STREQ(s1, b);
  EXPECT_EQ(strlen(s1), b.len());
}

class IOTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    tf = tmpfile();
    if(!tf)
      throw std::runtime_error("Can't create tmp file");
    l1 = "Hello there\n";
    l2 = "A very long and meaningless sentence-line.";
    fputs(l1, tf); 
    fputs(l2, tf);
    rewind(tf);
  }

  FILE *tf;
  const char *l1;
  const char *l2;
};

TEST_F(IOTest, fgets) {
  charb b(4);
  char *res = fgets(b, 3, tf); // 3 ignored!
  EXPECT_STREQ(l1, res);
  EXPECT_STREQ(l1, b);
  EXPECT_EQ(20, b.capacity());
  EXPECT_EQ(strlen(l1), b.len());

  res = fgets(b, tf);
  EXPECT_STREQ(l2, res);
  EXPECT_STREQ(l2, b);
  EXPECT_EQ(80, b.capacity());
  EXPECT_EQ(strlen(l2), b.len());

  res = fgets(b, tf);
  EXPECT_EQ((char *)0, res);
  EXPECT_STREQ(l2, b);
}

TEST_F(IOTest, getline) {
  charb b(5);
  ssize_t res;

  res = getline(b, tf);
  EXPECT_EQ(strlen(l1), res);
  EXPECT_STREQ(l1, b);
  EXPECT_LE(strlen(l1), b.len());

  res = getline(b, tf);
  EXPECT_EQ(strlen(l2), res);
  EXPECT_STREQ(l2, b);
  EXPECT_LE(strlen(l2), b.len());
  
  res = getline(b, tf);
  EXPECT_EQ(-1, res);
}

TEST_F(IOTest, fgets_append) {
  charb b(5);
  char *res;

  res = fgets_append(b, tf);
  EXPECT_STREQ(l1, res);
  EXPECT_STREQ(l1, b);

  res = fgets_append(b, tf);
  EXPECT_STREQ(l2, res);
  std::string str_res(l1);
  str_res += l2;
  EXPECT_STREQ(str_res.c_str(), b);
  EXPECT_EQ(strlen(l1) + strlen(l2), b.len());
}

TEST(CharbBasic, sprintf) {
  charb b(10);
  const char *fmt = "Hello %d times";
  const char *str_res = "Hello 1000 times";
  int res = sprintf(b, fmt, 1000);
  EXPECT_EQ(strlen(str_res), res);
  EXPECT_STREQ(str_res, b);
  EXPECT_EQ(strlen(str_res), b.len());
}

