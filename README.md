KMatLib
=======

Tiny matrix calculation library.

The features are as follows.

1. This library uses the stack instaed of the heap. 
2. This is a simple C++ template library. (in other words, This library is under development ).
3. In this library, a part written by me is public domain license. (Notice: Obviously Google C++ Unit Testing framework part is the new BSD License. )

By using the stack, high performace can gain compared with ordinary matrix library.
This library may adopt for a tiny matrix ( about N <= 10 ) because "stack overflow" may be occur. The stack is so limited resource, then this library cannot handle big matrix.
 
日本語での説明
=======
このライブラリは小さな行列に対する計算ライブラリです．

次の特徴を持ちます

1. スタックを用います．(ヒープではなく）
2. C++テンプレートライブラリでシンプルです．（言い換えれば，機能が限定されていて，開発途上です)
3. Public Domainライセンスです(もちろん私が書いた部分だけ)．利便性のためにGoogle C++ Testing Frameworkのソースを含めていますが，それは当然別ライセンスなのでご注意ください)

スタックを使うことで，通常の行列ライブラリよりは高い性能が得られます．しかし，スタックはヒープと比べるとかなり限られた資源です．気をつけないと，スタックオーバーフローが起こる可能性があります．このライブラリは，N<=10くらいの小さな行列に対してのみ適用することができるでしょう．

oniprog  tkuman@gmail.com

