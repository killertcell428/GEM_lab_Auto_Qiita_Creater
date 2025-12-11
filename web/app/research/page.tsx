'use client';

import { useEffect, useState } from 'react';
import { api } from '@/lib/api';
import LoadingSpinner from '@/components/LoadingSpinner';
import Link from 'next/link';
import { format } from 'date-fns';
import { ja } from 'date-fns/locale/ja';

interface ResearchItem {
  article_id: string;
  title: string;
  topic?: string;
  research_report?: string;
  has_research: boolean;
  has_plan: boolean;
  created_at?: string;
  updated_at?: string;
}

export default function ResearchLibraryPage() {
  const [researchItems, setResearchItems] = useState<ResearchItem[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [searchQuery, setSearchQuery] = useState('');

  useEffect(() => {
    loadResearchLibrary();
  }, []);

  const loadResearchLibrary = async () => {
    try {
      setLoading(true);
      setError(null);
      const data = await api.getResearchLibrary();
      setResearchItems(data);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'リサーチライブラリの読み込みに失敗しました';
      setError(errorMessage);
      console.error('リサーチライブラリ読み込みエラー:', err);
    } finally {
      setLoading(false);
    }
  };

  const filteredItems = researchItems.filter(item => {
    if (!searchQuery) return true;
    const query = searchQuery.toLowerCase();
    return (
      item.title.toLowerCase().includes(query) ||
      item.topic?.toLowerCase().includes(query) ||
      item.research_report?.toLowerCase().includes(query)
    );
  });

  return (
    <div className="space-y-6">
      <div className="flex justify-between items-center">
        <h1 className="text-3xl font-bold text-gray-900 dark:text-white">リサーチライブラリ</h1>
        <Link
          href="/"
          className="px-4 py-2 bg-gray-600 text-white rounded-md hover:bg-gray-700 transition-colors"
        >
          ダッシュボードに戻る
        </Link>
      </div>

      {error && (
        <div className="mb-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          {error}
        </div>
      )}

      {/* 検索バー */}
      <div>
        <input
          type="text"
          placeholder="タイトル、トピック、リサーチ結果で検索..."
          value={searchQuery}
          onChange={(e) => setSearchQuery(e.target.value)}
          className="w-full px-4 py-2 border border-gray-300 rounded-md dark:bg-gray-700 dark:border-gray-600 dark:text-white"
        />
      </div>

      {loading ? (
        <LoadingSpinner />
      ) : filteredItems.length === 0 ? (
        <div className="text-center py-12">
          <p className="text-gray-600 dark:text-gray-400 mb-4">
            {searchQuery ? '検索結果がありません' : 'リサーチ結果がありません'}
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
          {filteredItems.map((item) => (
            <Link
              key={item.article_id}
              href={`/research/${item.article_id}`}
              className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 hover:shadow-lg transition-shadow"
            >
              <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2 line-clamp-2">
                {item.title}
              </h3>
              {item.topic && (
                <p className="text-sm text-gray-600 dark:text-gray-400 mb-3 line-clamp-2">
                  {item.topic}
                </p>
              )}
              <div className="flex items-center space-x-2 mb-3">
                {item.has_research && (
                  <span className="px-2 py-1 text-xs font-medium bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200 rounded">
                    リサーチあり
                  </span>
                )}
                {item.has_plan && (
                  <span className="px-2 py-1 text-xs font-medium bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200 rounded">
                    プランあり
                  </span>
                )}
              </div>
              {item.research_report && (
                <p className="text-xs text-gray-500 dark:text-gray-500 line-clamp-3 mb-3">
                  {item.research_report.substring(0, 150)}...
                </p>
              )}
              <div className="text-xs text-gray-500 dark:text-gray-500">
                {item.created_at ? format(new Date(item.created_at), 'yyyy年MM月dd日', { locale: ja }) : '-'}
              </div>
            </Link>
          ))}
        </div>
      )}
    </div>
  );
}

