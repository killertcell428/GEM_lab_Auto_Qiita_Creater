'use client';

import { useEffect, useState } from 'react';
import { api, ArticleListItem } from '@/lib/api';
import ArticleStatusCard from '@/components/ArticleStatusCard';
import LoadingSpinner from '@/components/LoadingSpinner';
import Link from 'next/link';

export default function Dashboard() {
  const [articles, setArticles] = useState<ArticleListItem[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    loadArticles();
  }, []);

  const loadArticles = async () => {
    try {
      setLoading(true);
      setError(null);
      const data = await api.listArticles();
      setArticles(data);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : '記事の読み込みに失敗しました';
      setError(errorMessage);
      console.error('記事読み込みエラー:', err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div>
      <div className="flex justify-between items-center mb-8">
        <h1 className="text-3xl font-bold text-gray-900 dark:text-white">ダッシュボード</h1>
        <Link
          href="/new"
          className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
        >
          新規記事作成
        </Link>
      </div>

      {error && (
        <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          {error}
        </div>
      )}

      {loading ? (
        <LoadingSpinner />
      ) : articles.length === 0 ? (
        <div className="text-center py-12">
          <p className="text-gray-600 dark:text-gray-400 mb-4">記事がありません</p>
          <Link
            href="/new"
            className="inline-block px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
          >
            最初の記事を作成
          </Link>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {articles.map((article) => (
            <ArticleStatusCard key={article.id} article={article} />
          ))}
        </div>
      )}
    </div>
  );
}
